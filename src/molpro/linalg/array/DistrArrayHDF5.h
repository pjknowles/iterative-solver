#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DISTRARRAYHDF5_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DISTRARRAYHDF5_H
#include <hdf5.h>
#include <string>

#include "molpro/linalg/array/DistrArrayDisk.h"

namespace molpro::linalg::array {
namespace util {
class PHDF5Handle;
}

/*!
 * @brief Distributed array storing the buffer on disk using hdf5.
 *
 * @warning Before using operations that require access to the file, the hdf5 file must be opened.
 * To avoid corruption on accidental termination the file must be closed when the data is no longer needed.
 *
 * Use cases
 * ---------
 * Backing up:
 * @code{.cpp}
 * auto backup(const DistrArray &current_solution, std::shared_ptr<PHDF5Handle> backup_file){
 *      return DistrArrayDisk{current_solution, backup_file};
 * }
 * @endcode
 */
class DistrArrayHDF5 : public DistrArrayDisk {
protected:
  std::shared_ptr<util::PHDF5Handle> m_file_handle; //!< hdf5 file handle
  hid_t m_dataset = dataset_default;                //!< HDF5 dataset object
public:
  static constexpr int dataset_default = -1; //!< default value for dataset id
  const std::string dataset_name = "array";  //!< name of HDF5 dataset where array is stored

  DistrArrayHDF5() = delete;
  DistrArrayHDF5(const DistrArrayHDF5 &source) = delete;
  //! Takes ownership of source content.
  DistrArrayHDF5(DistrArrayHDF5 &&source) noexcept;

  /*!
   * @brief Create a disk array with a file assigned and default distribution.
   *
   * By default array is distributed uniformly with any remainder spread over the first processes so their size is
   * larger by 1.
   *
   * @note if the file contains dataset with dataset_name than it should be of the same dimension
   *
   * @param file_handle handle for opening the HDF5 group where array is/will be stored.
   * @param dimension size of array. If dataset already exists it will be resized to dimension.
   */
  DistrArrayHDF5(const std::shared_ptr<util::PHDF5Handle> &file_handle, size_t dimension);
  /*!
   * @brief Create a disk array with file and distribution assigned.
   *
   * @note if the file contains dataset with dataset_name than it should be of the same dimension
   *
   * @param file_handle handle for opening the HDF5 group where array is/will be stored.
   * @param distribution specifies how array is distributed among processes
   */
  DistrArrayHDF5(const std::shared_ptr<util::PHDF5Handle> &file_handle, std::unique_ptr<Distribution> distribution);
  /*!
   * @brief Creates a disk array by copying source to disk.
   * @param source a distributed array
   * @param file_handle handle for opening the HDF5 group where array is/will be stored.
   */
  DistrArrayHDF5(const DistrArray &source, const std::shared_ptr<util::PHDF5Handle> &file_handle);

  /*!
   * @brief Create a copy of source array using a temporary file which will be erased on destruction
   * @param source array to be copied
   * @param base_name base name for the temporary file. It will form the first part of the temporary array name. The
   * middle will be chosen to be unique and the suffix will be ".hdf5"
   */
  static DistrArrayHDF5 CreateTempCopy(const DistrArray &source, const std::string &base_name = ".temp_array");

  /*!
   * @brief Copies the array from source.
   *
   * If this array is a dummy (i.e. has no underlying array) than
   *
   * Copies to the memory view if it was allocated, otherwise copies directly from disk
   * if source has opened access.
   * If no file handle was assigned, than it is copied from source.
   * If source is not compatible, than this will be a direct copy of source.
   *
   * @todo This is very unintuitive. I am deleting this operator until the logic can be made clear. The basic case of
   * copying data is already covered by copy().
   *
   * @param source
   * @return
   */
  DistrArrayHDF5 &operator=(const DistrArrayHDF5 &source) = delete;
  DistrArrayHDF5 &operator=(DistrArrayHDF5 &&source) noexcept;

  friend void swap(DistrArrayHDF5 &x, DistrArrayHDF5 &y) noexcept;

  ~DistrArrayHDF5() override;

  bool compatible(const DistrArrayHDF5 &source) const;

  //! Removes link to the array dataset from the hdf5 file. This does not reduce the file size, consider using h5repack.
  void erase() override;
  value_type at(index_type ind) const override;
  void set(index_type ind, value_type val) override;
  void get(index_type lo, index_type hi, value_type *buf) const override;
  std::vector<value_type> get(index_type lo, index_type hi) const override;
  void put(index_type lo, index_type hi, const value_type *data) override;
  //! Accumulates data into array. @note There is not special accumulate functionality, implementation does a get() and
  //! a put(). This might be done more efficiently manually.
  void acc(index_type lo, index_type hi, const value_type *data) override;
  std::vector<value_type> gather(const std::vector<index_type> &indices) const override;
  void scatter(const std::vector<index_type> &indices, const std::vector<value_type> &data) override;
  void scatter_acc(std::vector<index_type> &indices, const std::vector<value_type> &data) override;
  std::vector<value_type> vec() const override;

  //! @returns a copy of a file handle
  std::shared_ptr<util::PHDF5Handle> file_handle() const;

  //! @returns id of dataset object.
  hid_t dataset() const;
  //! Checks that the array dataset exists in the file. The group must be already open.
  int dataset_exists() const;
  //! True if dataset is currently open, implies that file and group is open as well.
  bool dataset_is_open() const;

private:
  void open_access();
  void close_access();
};

} // namespace molpro::linalg::array
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DISTRARRAYHDF5_H
