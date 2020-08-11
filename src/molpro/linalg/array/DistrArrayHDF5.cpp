#include "DistrArrayHDF5.h"
#include <algorithm>
#include <cassert>

#include "PHDF5Handle.h"
namespace molpro {
namespace linalg {
namespace array {

int DistrArrayHDF5::dataset_exists() const {
  assert(m_file_handle->group_is_open());
  return util::hdf5_link_exists(m_file_handle->group_id(), dataset_name);
}

void DistrArrayHDF5::open_access() {
  assert(m_file_handle);
  m_file_handle->open_group();
  auto exists = dataset_exists();
  if (exists > 0) {
    m_dataset = H5Dopen(m_file_handle->group_id(), dataset_name.c_str(), H5P_DEFAULT);
  } else if (exists == 0) {
    hsize_t dimensions[1] = {m_dimension};
    auto space = H5Screate_simple(1, dimensions, nullptr);
    m_dataset = H5Dcreate(m_file_handle->group_id(), dataset_name.c_str(), H5T_NATIVE_DOUBLE, space, H5P_DEFAULT,
                          H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space);
  } else
    MPI_Abort(communicator(), 1); // dataset failed to open
}

void DistrArrayHDF5::close_access() {
  if (m_dataset != dataset_default)
    H5Dclose(m_dataset);
  m_dataset = dataset_default;
  m_file_handle->close_file();
}

void DistrArrayHDF5::erase() {
  open_access();
  auto lapl_id = H5Dget_access_plist(m_file_handle->group_id());
  H5Ldelete(m_file_handle->group_id(), dataset_name.c_str(), lapl_id);
  H5Pclose(lapl_id);
  close_access();
}

bool DistrArrayHDF5::dataset_is_open() const {
  return m_file_handle && m_file_handle->group_is_open() && m_dataset != dataset_default;
}

DistrArray::value_type DistrArrayHDF5::at(index_type ind) const {
  value_type val;
  get(ind, ind, &val);
  return val;
}

void DistrArrayHDF5::set(index_type ind, value_type val) { put(ind, ind, &val); }

void DistrArrayHDF5::get(index_type lo, index_type hi, value_type *buf) const {
  if (!dataset_is_open())
    error("call open_access() before RMA operations");
  hsize_t count[1] = {hi - lo + 1};
  hsize_t offset[1] = {lo};
  auto memspace = H5Screate_simple(1, count, nullptr);
  auto fspace = H5Dget_space(m_dataset);
  H5Sselect_hyperslab(fspace, H5S_SELECT_SET, offset, nullptr, count, nullptr);
  auto plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
  auto ierr = H5Dread(m_dataset, H5T_NATIVE_DOUBLE, memspace, fspace, plist_id, buf);
  if (ierr < 0)
    error("read failed");
}

std::vector<DistrArrayHDF5::value_type> DistrArrayHDF5::get(DistrArray::index_type lo,
                                                            DistrArray::index_type hi) const {
  auto buf = std::vector<value_type>(hi - lo + 1);
  get(lo, hi, &buf[0]);
  return buf;
}

void DistrArrayHDF5::put(DistrArray::index_type lo, DistrArray::index_type hi, const DistrArray::value_type *data) {
  if (!dataset_is_open())
    error("call open_access() before RMA operations");
  hsize_t count[1] = {hi - lo + 1};
  hsize_t offset[1] = {lo};
  auto memspace = H5Screate_simple(1, count, nullptr);
  auto fspace = H5Dget_space(m_dataset);
  H5Sselect_hyperslab(fspace, H5S_SELECT_SET, offset, nullptr, count, nullptr);
  auto plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
  auto ierr = H5Dwrite(m_dataset, H5T_NATIVE_DOUBLE, memspace, fspace, plist_id, data);
  if (ierr < 0)
    error("read failed");
}

void DistrArrayHDF5::acc(index_type lo, index_type hi, const value_type *data) {
  auto disk_copy = get(lo, hi);
  std::transform(begin(disk_copy), end(disk_copy), data, begin(disk_copy), [](auto &l, auto &r) { return l + r; });
  put(lo, hi, &disk_copy[0]);
}

std::vector<DistrArrayHDF5::value_type> DistrArrayHDF5::gather(const std::vector<index_type> &indices) const {
  if (!dataset_is_open())
    error("call open_access() before RMA operations");
  auto sz = indices.size();
  auto data = std::vector<value_type>(sz);
  hsize_t count[] = {sz};
  auto memspace = H5Screate_simple(1, count, nullptr);
  auto fspace = H5Dget_space(m_dataset);
  auto indices_copy = std::vector<hsize_t>(sz);
  std::copy(begin(indices), end(indices), begin(indices_copy));
  H5Sselect_elements(fspace, H5S_SELECT_SET, sz, &indices_copy[0]);
  auto plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
  auto ierr = H5Dread(m_dataset, H5T_NATIVE_DOUBLE, memspace, fspace, plist_id, &data[0]);
  if (ierr < 0)
    error("read failed");
  return data;
}

void DistrArrayHDF5::scatter(const std::vector<index_type> &indices, const std::vector<value_type> &data) {
  if (!dataset_is_open())
    error("call open_access() before RMA operations");
  auto sz = indices.size();
  hsize_t count[] = {sz};
  auto memspace = H5Screate_simple(1, count, nullptr);
  auto fspace = H5Dget_space(m_dataset);
  auto indices_copy = std::vector<hsize_t>(sz);
  std::copy(begin(indices), end(indices), begin(indices_copy));
  H5Sselect_elements(fspace, H5S_SELECT_SET, sz, &indices_copy[0]);
  auto plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
  auto ierr = H5Dwrite(m_dataset, H5T_NATIVE_DOUBLE, memspace, fspace, plist_id, &data[0]);
  if (ierr < 0)
    error("read failed");
}

void DistrArrayHDF5::scatter_acc(std::vector<index_type> &indices, const std::vector<value_type> &data) {
  auto disk_copy = gather(indices);
  std::transform(begin(data), end(data), begin(disk_copy), begin(disk_copy), [](auto &l, auto &r) { return l + r; });
  scatter(indices, disk_copy);
}

std::vector<DistrArrayHDF5::value_type> DistrArrayHDF5::vec() const { return get(0, m_dimension - 1); }

} // namespace array
} // namespace linalg
} // namespace molpro
