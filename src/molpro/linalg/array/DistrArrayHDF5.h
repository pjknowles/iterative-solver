#ifndef GCI_SRC_MOLPRO_GCI_ARRAY_DISTRARRAYHDF5_H
#define GCI_SRC_MOLPRO_GCI_ARRAY_DISTRARRAYHDF5_H
#include <hdf5.h>

#include "molpro/linalg/array/DistrArrayDisk.h"

namespace molpro {
namespace linalg {
namespace array {
namespace util {
class PHDF5Handle;
}

/*!
 * @brief Distributed array storing the buffer on disk using hdf5.
 *
 * @warning Before using operations that require access to the file, the hdf5 file must be opened.
 * To avoid corruption on accidental termination the file must be closed when the data is no longer needed.
 */
class DistrArrayHDF5 final : public DistrArrayDisk {
protected:
  std::shared_ptr<util::PHDF5Handle> m_file_handle; //!< hdf5 file handle
  hid_t m_dataset = dataset_default;                //!< HDF5 dataset object
  const hid_t dataset_default = -1;                 //!< default value for dataset id
public:
  const std::string dataset_name = "array"; //!< name of HDF5 dataset where array is stored

  //! Constructor for a blank object. The blank is only useful as a temporary. Move a valid object inside the blank to
  //! make it usable.
  DistrArrayHDF5() = default;
  //! Copies the array from source. Copies to the memory view if it was allocated, otherwise copies directly from disk.
  DistrArrayHDF5(const DistrArrayHDF5 &source);
  //! Takes ownership of source content.
  DistrArrayHDF5(DistrArrayHDF5 &&source);

  /*!
   * @brief
   * @param file_handle handle for opening the HDF5 group where array is/will be stored.
   * @param dimension size of array. If dataset already exists it will be resized to dimension.
   * @param prof profiler
   */
  DistrArrayHDF5(std::shared_ptr<util::PHDF5Handle> file_handle, size_t dimension,
                 std::shared_ptr<Profiler> prof = nullptr);
  /*!
   * @brief
   * @param file_handle
   * @param prof
   */
  DistrArrayHDF5(std::shared_ptr<util::PHDF5Handle> file_handle, std::shared_ptr<Profiler> prof = nullptr);

  DistrArrayHDF5 &operator=(const DistrArrayHDF5 &);
  DistrArrayHDF5 &operator=(DistrArrayHDF5 &&);

  DistrArrayHDF5 &operator=(const DistrArray &source);

  friend void swap(DistrArrayHDF5 &x, DistrArrayHDF5 &y);

  ~DistrArrayHDF5() override;

  //! @returns HDF5 dataset object
  hid_t datest() const { return m_dataset; }

protected:
  class LocalBufferDiskHDF5 : public DistrArrayHDF5::LocalBufferDisk {
  public:
    explicit LocalBufferDiskHDF5(DistrArrayHDF5 &source);
  };

public:
  [[nodiscard]] const Distribution &distribution() const override;

  void open_access() override;
  void close_access() override;
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

protected:
  //! Checks that the array dataset under default name exists. The group must be already open.
  int dataset_exists() const;
  bool dataset_is_open() const;
};

} // namespace array
} // namespace linalg
} // namespace molpro
#endif // GCI_SRC_MOLPRO_GCI_ARRAY_DISTRARRAYHDF5_H
