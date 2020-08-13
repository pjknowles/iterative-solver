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

class DistrFlags;
} // namespace util
/*!
 * @brief Distributed array located primarily on disk
 *
 * This class stores the full array on disk. RMA operations read data from disk. A memory view buffer can be allocated
 * which will read the local buffer from disk and keep array in memory until it is flushed. This allows multiple
 * operations on local section to be performed without I/O.
 *
 * Concepts
 * --------
 *   * memory view - a buffer in RAM which mirrors the content of local section of the array on disk.
 *   * a snapshot - a temporary local buffer with contents read directly from disk.
 *
 * Differences with DistrArray
 * ---------------------------
 *
 * open_access() opens access to the underlying data on disk and close_access() closes that access. If the underlying
 * storage space does not exist than open_access() creates it and reserves necessary space (e.g. creates a file or
 * hdf5 dataset). Both functions are assumed to be collective, although not all implementations might require that.
 *
 * allocate_buffer() creates the memory view buffer. The buffer will persist until free_buffer() is called. The memory
 * view has to be written to disk by calling flush.
 *
 * LocalBuffer will use the memory view buffer if it was allocated, otherwise it will allocate a new buffer. The new
 * buffer will be a snapshot of the array currently on disk. A non-const local buffer snapshot is always written to disk
 * on destruction, even if it was not modified.
 *
 * Linear algebra operations work by creating a LocalBuffer instance which loads relevant section of the array,
 * followed by numerical operations and optionally dumping of LocalBuffer to disk.
 *
 * There is no default mechanism for keeping memory view in sync. That is, modifications of array on disk via RMA calls
 * makes memory view out of sync.
 *
 * @note The main point of memory view is to minimize IO for linear algebra operations.
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
 * @warning Mixing threaded read and write operations can lead to undefined behaviour
 * @warning DistrArrayDisk is **NOT** thread safe by itself, even though it uses threads for I/O. That is, multiple
 * threads calling DistrArrayDisk routines can lead to undefined behavior.
 *
 * @code{.cpp}
 * auto da = DistrArrayDisk{...};
 * da.open_access();
 * ...
 * da.put(lo,hi, data); // Normal put as in base class is guaranteed to finish
 * da.tput(lo, hi, data); // same as put, because DiskIO syncs on destruction
 * auto t = da.tput(); // does I/O in a thread
 * // do something time consuming
 * t.wait(); // wait for the thread to finish the I/O
 * @endcode
 *
 * Advice
 * ------
 *
 * Enclose section of linear algebra calls with  allocate_buffer(), flush() and free_buffer().
 * Operations will be done using memory view buffer without unnecessary IO.
 *
 * Mixing up incompatible RMA and linear algebra operations can lead to undefined behavior of DistrArray. This is also
 * true for DistrArrayDisk and the buffering mechanism where the allocated buffer mimics data stored on disk.
 *
 * @code{.cpp}
 * auto x = DistrArrayDisk{...};
 * auto y = DistrArrayDisk{...};
 * x.open_access();
 * y.open_access();
 * x.allocate_buffer();
 * y.allocate_buffer();
 * auto norm = std::sqrt(x.dot(x));
 * auto ov = x.dot(y);
 * y.axpy(ov / norm, x);
 * y.flush();
 * y.free_buffer();
 * x.free_buffer();
 * x.close_access();
 * y.close_access();
 * @endcode
 *
 *
 */
class DistrArrayDisk : public DistrArray {
protected:
  bool m_allocated = false;       //!< Flags that the memory view buffer has been allocated
  Span<value_type> m_view_buffer; //!< memory view buffer either wraps allocated buffer or stores user supplied buffer
  std::vector<value_type> m_owned_buffer; //!< buffer allocated by the class
  using DistrArray::DistrArray;

  DistrArrayDisk();
  DistrArrayDisk(const DistrArrayDisk &source);
  DistrArrayDisk(DistrArrayDisk &&source) noexcept;

public:
  //! Opens access to the storage on disk, creating the underlying storage object if it does not exist. Assume
  //! collective, but not all implementations might require that.
  virtual void open_access() = 0;
  //! Close access to the storage on disk, buffer is flushed to disk. Assume collective, but not all implementations
  //! might require that.
  virtual void close_access() = 0;
  void allocate_buffer() override;
  /*!
   * @brief  Uses a provided buffer.
   *
   * If a buffer has already been allocated than its content is copied into the new buffer. Any old LocalBuffer
   * instances will still be using the old buffer.
   *
   * @param buffer memory view buffer to use.
   */
  void allocate_buffer(Span<value_type> buffer);
  bool empty() const override;
  //! Release the allocated buffer. @note buffer is not flushed.
  void free_buffer() override;
  //! Writes the memory view buffer.
  virtual void flush();
  //! Erase the array from disk.
  virtual void erase() = 0;

protected:
  class LocalBufferDisk : public DistrArray::LocalBuffer {
  public:
    explicit LocalBufferDisk(DistrArrayDisk &source);
    ~LocalBufferDisk() override;
    //! Prints true if the buffer is a snapshot and not a memory view
    bool is_snapshot();
    //! access dumping flag. If true, than buffer is dumped to file on destruction.
    bool &dump() { return m_dump; };
    const bool &dump() const { return m_dump; };

  protected:
    bool m_dump = true;                        //! If true than snapshot buffer is dumped to disk on destruction
    std::vector<value_type> m_snapshot_buffer; //!< buffer for storing a snap shot when memory view is not allocated
    DistrArrayDisk &m_source;                  //!< keep a handle on source to dump data to disk
  };

public:
  [[nodiscard]] std::unique_ptr<LocalBuffer> local_buffer() override;
  [[nodiscard]] std::unique_ptr<const LocalBuffer> local_buffer() const override;

  //! Writes section of array to disk in a separate thread, bypassing the memory view.
  //! @return Task manager that can test and/or wait for completion
  std::unique_ptr<util::Task<void>> tput(index_type lo, index_type hi, const value_type *data);
  /*!
   * @brief Reads section of array from disk in a separate thread
   * @warning Reading and writing bypasses the memory-view, the buffer should be flushed before
   * @return Task manager that can test and/or wait for completion
   */
  std::unique_ptr<util::Task<void>> tget(index_type lo, index_type hi, value_type *buf);
  std::unique_ptr<util::Task<std::vector<value_type>>> tget(index_type lo, index_type hi);
  std::unique_ptr<util::Task<value_type>> tat(index_type ind) const;
  std::unique_ptr<util::Task<void>> tset(index_type ind, value_type val);
  std::unique_ptr<util::Task<void>> tacc(index_type lo, index_type hi, const value_type *data);
  std::unique_ptr<util::Task<std::vector<value_type>>> tgather(const std::vector<index_type> &indices) const;
  std::unique_ptr<util::Task<void>> tscatter(const std::vector<index_type> &indices,
                                             const std::vector<value_type> &data);
  std::unique_ptr<util::Task<void>> tscatter_acc(std::vector<index_type> &indices, const std::vector<value_type> &data);
  std::unique_ptr<util::Task<std::vector<value_type>>> tvec() const;

  //! Loads the local buffer in a separate thread
  std::unique_ptr<util::Task<std::unique_ptr<LocalBuffer>>> tlocal_buffer();
  std::unique_ptr<util::Task<std::unique_ptr<const LocalBuffer>>> tlocal_buffer() const;
};

} // namespace array
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DISTRARRAYDISK_H
