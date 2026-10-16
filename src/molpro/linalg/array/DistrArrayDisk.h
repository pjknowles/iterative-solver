#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DISTRARRAYDISK_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DISTRARRAYDISK_H

#include <future>
#include <molpro/Profiler.h>
#include <molpro/linalg/array/DistrArray.h>
#include <molpro/linalg/array/Span.h>
#include <string>

namespace molpro::linalg::array {
/*!
 * @brief Distributed array located primarily on disk
 *
 * This class stores the full array on disk and implements RMA and more efficient linear algebra operations.
 *
 * RMA operations read/write directly to disk.
 *
 * There is no local buffer: the memory used must not scale with the size of the array, so local_buffer() throws.
 * Elementwise operations, including those with arrays held in memory, page through the local section
 * disk_page_size() elements at a time.
 *
 * BufferManager reads the local section in chunks using a separate thread for I/O. This is more memory efficient
 * and allows overlap of communication and computation. The walk through the local section is via iterators.
 *
 * IO can be done in a separate thread using util::Task.
 *
 * @code{.cpp}
 * #include <molpro/linalg/array/util.h>
 * auto da = DistrArrayDisk{...};
 * da.put(lo,hi, data); // Normal put same as in the base class is guaranteed to finish
 * auto t = util::Task::create([&](){da.put(lo, hi, data)}); // does I/O in a new thread
 * // do something time consuming
 * t.wait(); // wait for the thread to finish the I/O
 * @endcode
 *
 */
class DistrArrayDisk : public DistrArray {
public:
  using disk_array = void; //!< a compile time tag that this is a distributed disk array
protected:
  bool m_allocated = false;                     //!< Flags that the memory view buffer has been allocated
  std::unique_ptr<Distribution> m_distribution; //!< describes distribution of array among processes
  size_t m_buffer_size = 8192;                  //!< buffer size for paged access via BufferManager
  using DistrArray::DistrArray;

  DistrArrayDisk(std::unique_ptr<Distribution> distr, MPI_Comm commun);
  DistrArrayDisk();
  DistrArrayDisk(const DistrArray& source);
  DistrArrayDisk(DistrArrayDisk&& source) noexcept;
  ~DistrArrayDisk() override;

public:
  //! Erase the array from disk.
  virtual void erase() = 0;
  [[nodiscard]] const Distribution& distribution() const override;
  [[nodiscard]] value_type dot(const DistrArrayDisk& y) const { return DistrArray::dot(y); }
  using DistrArray::dot;
  void set_buffer_size(size_t buffer_size) { m_buffer_size = buffer_size; }
  [[nodiscard]] size_t disk_page_size() const override { return m_buffer_size; }

  //! The directory of the file holding this array, or an empty string if that is not known
  [[nodiscard]] virtual std::string file_directory() const { return {}; }

  /*!
   * @brief The directory for the file of a temporary copy of source
   * @param source the array to be copied
   * @param directory the directory asked for, or an empty string if none was
   * @return directory if not empty, otherwise the directory of the file of source if it is a disk array, otherwise
   * the current directory
   */
  static std::string temp_copy_directory(const DistrArray& source, const std::string& directory);

  //! Not available: an array on disk is not held in memory. @throws std::logic_error
  [[nodiscard]] std::unique_ptr<LocalBuffer> local_buffer() override;
  //! Not available: an array on disk is not held in memory. @throws std::logic_error
  [[nodiscard]] std::unique_ptr<const LocalBuffer> local_buffer() const override;
};

double dot(const DistrArrayDisk& x, const DistrArrayDisk& y);
double dot(const DistrArrayDisk& x, const DistrArray& y);
double dot(const DistrArray& x, const DistrArrayDisk& y);

} // namespace molpro::linalg::array

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DISTRARRAYDISK_H
