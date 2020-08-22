#include "DistrArrayHDF5.h"
#include "util/Distribution.h"
#include <algorithm>

#include "PHDF5Handle.h"
namespace molpro {
namespace linalg {
namespace array {
namespace {
int mpi_size(MPI_Comm comm) {
  int rank;
  MPI_Comm_size(comm, &rank);
  return rank;
}
} // namespace

DistrArrayHDF5::DistrArrayHDF5(std::shared_ptr<util::PHDF5Handle> file_handle,
                               std::unique_ptr<Distribution> distribution, std::shared_ptr<Profiler> prof)
    : DistrArrayDisk(std::move(distribution), file_handle->communicator(), std::move(prof)),
      m_file_handle(std::move(file_handle)) {
  if (m_distribution->border().first != 0)
    DistrArray::error("Distribution of array must start from 0");
}

DistrArrayHDF5::DistrArrayHDF5(std::shared_ptr<util::PHDF5Handle> file_handle, size_t dimension,
                               std::shared_ptr<Profiler> prof)
    : DistrArrayHDF5(std::move(file_handle),
                     std::make_unique<Distribution>(util::make_distribution_spread_remainder<index_type>(
                         dimension, mpi_size(file_handle->communicator()))),
                     std::move(prof)) {}

DistrArrayHDF5::DistrArrayHDF5(std::shared_ptr<util::PHDF5Handle> file_handle, std::shared_ptr<Profiler> prof)
    : DistrArrayDisk(std::make_unique<Distribution>(std::vector<index_type>{0, 0}), file_handle->communicator(),
                     std::move(prof)),
      m_file_handle(std::move(file_handle)) {
  m_file_handle->open_group();
  if (m_file_handle->group_is_open() && dataset_exists()) {
    DistrArrayHDF5::open_access();
    auto space = H5Dget_space(m_dataset);
    auto n = H5Sget_simple_extent_ndims(space);
    if (n < 0)
      DistrArray::error("failed to get number of dimensions in dataset");
    if (n != 1)
      DistrArray::error("dataset does not correspond to a 1D array");
    hsize_t dims;
    n = H5Sget_simple_extent_dims(space, &dims, nullptr);
    if (n < 0)
      DistrArray::error("failed to get dimensions in dataset");
    // create a new array
    DistrArrayHDF5 t{m_file_handle, dims, std::move(prof)};
    swap(*this, t);
    t.m_file_handle.reset();
  }
  m_file_handle->close_file();
}

DistrArrayHDF5::DistrArrayHDF5(const DistrArray &source, std::shared_ptr<util::PHDF5Handle> file_handle)
    : DistrArrayDisk(std::make_unique<Distribution>(source.distribution()), source.communicator(), source.m_prof),
      m_file_handle(std::move(file_handle)) {
  if (!source.empty()) {
    DistrArrayHDF5::open_access();
    DistrArrayHDF5::copy(source);
    DistrArrayHDF5::close_access();
  }
}

DistrArrayHDF5::DistrArrayHDF5() = default;

DistrArrayHDF5::DistrArrayHDF5(const DistrArrayHDF5 &source)
    : DistrArrayDisk(source), m_file_handle(source.m_file_handle), m_dataset(source.m_dataset) {}

DistrArrayHDF5::DistrArrayHDF5(DistrArrayHDF5 &&source) noexcept
    : DistrArrayDisk(std::move(source)), m_file_handle(std::move(source.m_file_handle)), m_dataset(source.m_dataset) {
  source.m_dataset = source.dataset_default;
}

DistrArrayHDF5 &DistrArrayHDF5::operator=(DistrArrayHDF5 &&source) noexcept {
  DistrArrayHDF5 t{std::move(source)};
  swap(*this, t);
  return *this;
}

DistrArrayHDF5::~DistrArrayHDF5() {
  if (m_allocated && dataset_is_open())
    DistrArrayHDF5::flush();
  if (dataset_is_open())
    DistrArrayHDF5::close_access();
}

bool DistrArrayHDF5::compatible(const DistrArrayHDF5 &source) const {
  auto res = DistrArray::compatible(source);
  if (m_distribution && source.m_distribution)
    res &= m_distribution->compatible(*source.m_distribution);
  else
    res &= !m_distribution && !source.m_distribution;
  return res;
}

void swap(DistrArrayHDF5 &x, DistrArrayHDF5 &y) noexcept {
  using std::swap;
  swap(x.m_dimension, y.m_dimension);
  swap(x.m_communicator, y.m_communicator);
  swap(x.m_prof, y.m_prof);
  swap(x.m_allocated, y.m_allocated);
  swap(x.m_view_buffer, y.m_view_buffer);
  swap(x.m_owned_buffer, y.m_owned_buffer);
  swap(x.m_distribution, y.m_distribution);
  swap(x.m_file_handle, y.m_file_handle);
  swap(x.m_dataset, y.m_dataset);
}

int DistrArrayHDF5::dataset_exists() const {
  if (!m_file_handle->group_is_open())
    error("group must be open before using dataset");
  return util::hdf5_link_exists(m_file_handle->group_id(), dataset_name);
}

void DistrArrayHDF5::open_access() {
  if (dataset_is_open())
    return;
  if (!m_file_handle)
    error("must provide a file handle before openning access to disk array");
  auto id = m_file_handle->open_file(util::HDF5Handle::Access::read_write);
  if (id < 0)
    error("Failed to open file");
  id = m_file_handle->open_group();
  if (id < 0)
    error("Failed to open group");
  auto exists = dataset_exists();
  if (exists > 0) {
    if (m_dataset == dataset_default)
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

bool DistrArrayHDF5::empty() const {
  auto res = DistrArrayDisk::empty() && !dataset_is_open();
  return res;
}

void DistrArrayHDF5::erase() {
  open_access();
  auto lapl_id = H5Dget_access_plist(m_dataset);
  if (lapl_id < 0)
    error("Failed to access property list of a group");
  H5Ldelete(m_file_handle->group_id(), dataset_name.c_str(), lapl_id);
  H5Pclose(lapl_id);
  close_access();
}

bool DistrArrayHDF5::dataset_is_open() const {
  return m_file_handle && m_file_handle->group_is_open() && m_dataset != dataset_default;
}

DistrArray::value_type DistrArrayHDF5::at(index_type ind) const {
  value_type val;
  get(ind, ind + 1, &val);
  return val;
}

void DistrArrayHDF5::set(index_type ind, value_type val) { put(ind, ind + 1, &val); }

void DistrArrayHDF5::get(index_type lo, index_type hi, value_type *buf) const {
  if (lo >= hi)
    return;
  if (!dataset_is_open())
    error("call open_access() before RMA operations");
  hsize_t count[1] = {hi - lo};
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
  if (lo >= hi)
    return {};
  auto buf = std::vector<value_type>(hi - lo);
  get(lo, hi, &buf[0]);
  return buf;
}

void DistrArrayHDF5::put(DistrArray::index_type lo, DistrArray::index_type hi, const DistrArray::value_type *data) {
  if (!dataset_is_open())
    error("call open_access() before RMA operations");
  hsize_t count[1] = {hi - lo};
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

std::vector<DistrArrayHDF5::value_type> DistrArrayHDF5::vec() const { return get(0, m_dimension); }

std::shared_ptr<util::PHDF5Handle> DistrArrayHDF5::file_handle() const { return m_file_handle; }

hid_t DistrArrayHDF5::dataset() const { return m_dataset; }

} // namespace array
} // namespace linalg
} // namespace molpro
