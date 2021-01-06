#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DISTRARRAYDISK_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DISTRARRAYDISK_H

#include "molpro/linalg/array/DistrArray.h"
#include "molpro/linalg/array/Span.h"

namespace molpro::linalg::array {
namespace util {
template <typename Result>
class Task;

} // namespace util
/*!
 * @brief Distributed array located primarily on disk
 *
 * This class stores the full array on disk and implements RMA and more efficient linear algebra operations.
 *
 * RMA operations read/write directly to disk. Care must be taken that they do not overlap with local buffer
 * modifications.
 *
 * LocalBuffer reads the whole local section into memory. This might be prohivitively expensive for large arrays, so
 * care must be taken.
 *
 * LocalBufferChunked reads the local section in chunks using a separate thread for I/O. This is more memory efficient
 * and allows overlap of communication and computation.
 *
 * More efficient linear algebra operations are implemented using ChunkedLocalBuffer.
 *
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
 * @code{.cpp}
 * @endcode
 *
 *
 */
class DistrArrayDisk : public DistrArray {
public:
  using disk_array = void; //!< a compile time tag that this is a distributed disk array
protected:
  bool m_allocated = false;                     //!< Flags that the memory view buffer has been allocated
  std::unique_ptr<Distribution> m_distribution; //!< describes distribution of array among processes
  using DistrArray::DistrArray;

  DistrArrayDisk(std::unique_ptr<Distribution> distr, MPI_Comm commun);
  DistrArrayDisk();
  DistrArrayDisk(const DistrArrayDisk &source);
  DistrArrayDisk(DistrArrayDisk &&source) noexcept;
  ~DistrArrayDisk() override;

public:
  //! Erase the array from disk.
  virtual void erase() = 0;
  [[nodiscard]] const Distribution &distribution() const override;

protected:
  //! Reads the whole local buffer from disk into memory. By default the buffer is written to disk on destruction,
  //! unless do_dump is false.
  class LocalBufferDisk : public DistrArray::LocalBuffer {
  public:
    explicit LocalBufferDisk(DistrArrayDisk &source);
    explicit LocalBufferDisk(DistrArrayDisk &source, const span::Span<value_type> &buffer);
    ~LocalBufferDisk() override;
    //! If true, than buffer is dumped to file on destruction.
    bool do_dump = true;

  protected:
    std::vector<value_type> m_snapshot_buffer; //!< when external buffer is not provided
    DistrArrayDisk &m_source;                  //!< keep a handle on source to dump data to disk
  };

  // TODO The mechanism for loading the file in chunks should be shared with serial disk array
  //! Loads local buffer in small chunks
  class LocalBufferChunked {};

public:
  [[nodiscard]] std::unique_ptr<LocalBuffer> local_buffer() override;
  [[nodiscard]] std::unique_ptr<const LocalBuffer> local_buffer() const override;
  //! Access local section, reading it into the provided buffer
  [[nodiscard]] std::unique_ptr<LocalBuffer> local_buffer(const span::Span<value_type> &buffer);
  //! Read-only access to the local section, reading it into the provided buffer
  [[nodiscard]] std::unique_ptr<const LocalBuffer> local_buffer(const span::Span<value_type> &buffer) const;
  [[nodiscard]] std::unique_ptr<LocalBufferChunked> local_buffer_chunked();
  [[nodiscard]] std::unique_ptr<const LocalBufferChunked> local_buffer_chunked() const;

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

double dot(const DistrArrayDisk &x, const DistrArrayDisk &y);
double dot(const DistrArrayDisk &x, const DistrArray &y);
double dot(const DistrArray &x, const DistrArrayDisk &y);

} // namespace molpro::linalg::array

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DISTRARRAYDISK_H
