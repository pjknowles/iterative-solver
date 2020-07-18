#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DISTRARRAYDISK_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DISTRARRAYDISK_H

#include "molpro/linalg/array/DistrArray.h"
#include "molpro/linalg/array/Span.h"

namespace molpro {
namespace linalg {
namespace array {
namespace util {
template <typename Result>
class Task;
} // namespace util
/*!
 * @brief Distributed array located primarily on disk
 *
 * This class stores the full array on disk. RMA operations read data from disk. A memory view buffer can be allocated
 * which will read the local buffer from disk and keep array in memory until it is flushed. This allows multiple
 * operations on local section to be performed without I/O at the end of each one.
 *
 * Concepts
 * --------
 *   * memory view - a buffer in RAM which mirrors the content of local section of the array on disk.
 *   * a snapshot - a temporary buffer with contents read directly from disk, ignoring the memory view.
 *
 * Mechanical differences with DistrArray
 * --------------------------------------
 *
 * open_access() opens access to the underlying data on disk and close_access() closes that access. If the underlying
 * storage space does not exist than open_access() creates it and reserves necessary space (e.g. creates a file or
 * hdf5 dataset). Both functions are assumed to be collective, although not all implementations might require that.
 *
 * allocate_buffer() creates the memory view buffer. The buffer will persist until deallocate_buffer() is called.
 *
 * LocalBuffer will use the memory view buffer if it was allocated, otherwise it will allocate a new buffer. The new
 * buffer will be a snapshot of the array currently on disk. A non-const local buffer snapshot is always written to disk
 * on destruction, even if it was not modified.
 *
 * Linear algebra operations work by creating a LocalBuffer instance which loads relevant section of the array,
 * followed by numerical operations and optionally dumping of LocalBuffer to disk if it was modified.
 *
 * To minimize I/O the memory view buffer is reused as long as it is in sync with the disc version. Operations such as
 * writing to a non-local section of an array invalidates the local buffer on processes that own it. When the memory
 * view buffer is invalidated, it will be read from disk again on next usage.
 *
 * Flushing writes the memory view to disk if it was allocated and is valid.
 *
 * Synchronisation has a dual effect of synchronising the processes, and if the memory view was allocated synchronising
 * it with disk. If memory view is valid, than it is written to disk, otherwise it is read from disk (may be different
 * for each process).
 *
 * IO can be done with a thread
 * ----------------------------
 *
 * Some RMA operations have a threaded implementation their names are prefixed with "t", e.g tput(), tget(). The
 * underlying IO operation is done in a separate thread allowing user to overlap it with computation. The separate
 * thread is managed by util::Task.
 *
 * tlocal_buffer() also uses a thread for loading the data from disk. If the LocalBuffer uses the memory view than the
 * task will be completed immediately.
 *
 * @warning DistrArrayDisk is **NOT** thread safe by itself, even though it uses threads for I/O. That is, multiple
 * threads calling DistrArrayDisk routines can lead to undefined behavior.
 *
 * Advice
 * ------
 *
 * Call allocate_buffer() and deallocate_buffer() between a section of computation with many linear algebra calls.
 * Operations will be done using memory view buffer without any IO.  At the end a flush() or sync() call has to be made
 * to write the new array to disk.
 *
 * Mixing up incompatible RMA and linear algebra operations can lead to undefined behavior of DistrArray. This is also
 * true for DistrArrayDisk and the buffering mechanism where the allocated buffer mimics data stored on disk.
 *
 * @code{.cpp}
 * auto da = DistrArrayDisk{...};
 * ...
 * da.put(lo,hi, data); // Normal put as in base class is guaranteed to finish
 * da.tput(lo, hi, data); // same as put, because DiskIO syncs on destruction
 * auto t = da.tput(); // does I/O in a thread
 * // do something time consuming
 * t.wait(); // wait for the thread to finish the I/O
 * @endcode
 */
class DistrArrayDisk : public DistrArray {
protected:
  bool m_allocated;          //!< Flags that the memory view buffer has been allocated
  bool m_invalid;            //!< Flags that the allocated buffer is out of sync with the array on disk
  Span<value_type> m_buffer; //!< memory view buffer either wraps allocated buffer or stores user supplied buffer
  std::vector<value_type> m_allocated_buffer; //!< buffer allocated by the class
  explicit DistrArrayDisk(const DistrArrayDisk &source) = default;

public:
  //! Opens access to the storage on disk, creating the underlying storage object if it does not exist. Assume
  //! collective, but not all implementations might require that.
  virtual void open_access() = 0;
  //! Close access to the storage on disk, buffer is flushed to disk. Assume collective, but not all implementations
  //! might require that.
  virtual void close_access() = 0;
  void allocate_buffer() override;
  //! Uses provided buffer. @returns false if allocation fails because a buffer was already allocated
  bool allocate_buffer(Span<value_type> buffer);
  //! Release the allocated buffer
  virtual void deallocate_buffer();
  //! Writes memory view buffer
  virtual void flush();
  void sync() const override;
  bool compatible(const DistrArrayDisk &other) const;
  //! Erase the array from disk deleting the underlying storage object.
  virtual bool erase() = 0;
  //! Returns true if the allocated buffer is out of sync with the array on disk
  bool invalid() const { return m_invalid; };
  //! Marks memory view buffer out of sync
  void invalidate() { m_invalid = true; };

protected:
  class LocalBufferDisk : public DistrArray::LocalBuffer {
  public:
    LocalBufferDisk(DistrArrayDisk &source);

  protected:
    bool m_dump_on_destruct = true; //!< whether to dump contents of local buffer to disk on destruction
  };

public:
  [[nodiscard]] std::unique_ptr<LocalBuffer> local_buffer() override;
  [[nodiscard]] std::unique_ptr<const LocalBuffer> local_buffer() const override;
  [[nodiscard]] const Distribution &distribution() const override;

  value_type at(index_type ind) const override;
  void set(index_type ind, value_type val) override;
  std::vector<value_type> get(index_type lo, index_type hi) const override;
  void get(index_type lo, index_type hi, value_type *buf) const override;
  void put(index_type lo, index_type hi, const value_type *data) override;
  void acc(index_type lo, index_type hi, const value_type *data) override;
  std::vector<value_type> gather(const std::vector<index_type> &indices) const override;
  void scatter(const std::vector<index_type> &indices, const std::vector<value_type> &data) override;
  void scatter_acc(std::vector<index_type> &indices, const std::vector<value_type> &data) override;
  std::vector<value_type> vec() const override;

  //! Writes section of array to disk in a separate thread, bypassing the memory view.
  //! @return Task manager that can test and/or wait for completion
  virtual std::unique_ptr<util::Task<void>> tput() = 0;
  //!
  /*!
   * @brief Reads section of array from disk in a separate thread
   * @warning Reading and writing bypasses the memory-view, the buffer should be flushed before
   * @return Task manager that can test and/or wait for completion
   */
  virtual std::unique_ptr<util::Task<void>> tget() = 0;

protected:
public:
  //! Loads the local buffer in a separate thread
  std::pair<LocalBuffer, std::unique_ptr<util::Task<void>>> tlocal_buffer();
};

} // namespace array
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DISTRARRAYDISK_H
