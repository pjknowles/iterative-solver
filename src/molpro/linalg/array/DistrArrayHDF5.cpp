#include "DistrArrayHDF5.h"
#include "util/Distribution.h"
#include <algorithm>

#include "PHDF5Handle.h"
#include "util/TempHandle.h"
#include "util/fs.h"
#include "util/temp_phdf5_handle.h"
namespace molpro::linalg::array {
namespace {
int mpi_size(MPI_Comm comm) {
  int rank;
  MPI_Comm_size(comm, &rank);
  return rank;
}

//! A handle to a new temporary file, named after base_name, or in base_name if that is an existing directory
std::shared_ptr<util::PHDF5Handle> temp_file_handle(const std::string &base_name, MPI_Comm comm) {
  const auto name = fs::is_directory(base_name) ? (fs::path(base_name) / ".temp_array").string() : base_name;
  auto handle = std::make_shared<util::PHDF5Handle>(util::temp_phdf5_handle(name, comm));
  handle->assign_group("/");
  return handle;
}

//! Reports an error through DistrArray::error() if an HDF5 call returned a negative value, which signals failure
template <typename T>
T check(const DistrArray &array, T result, const std::string &operation) {
  if (result < 0)
    array.error("DistrArrayHDF5: failed to " + operation);
  return result;
}

//! Owns an HDF5 identifier and closes it on destruction, so that it is released when DistrArray::error() throws
class Hid {
public:
  Hid(const DistrArray &array, hid_t id, herr_t (*close)(hid_t), const std::string &operation)
      : m_id(check(array, id, operation)), m_close(close) {}
  Hid(Hid &&source) noexcept : m_id(source.m_id), m_close(source.m_close) { source.m_id = H5I_INVALID_HID; }
  Hid(const Hid &) = delete;
  Hid &operator=(const Hid &) = delete;
  Hid &operator=(Hid &&) = delete;
  ~Hid() {
    if (m_id >= 0)
      m_close(m_id);
  }
  operator hid_t() const { return m_id; }

private:
  hid_t m_id;
  herr_t (*m_close)(hid_t);
};

//! Transfer properties for independent parallel I/O
Hid independent_transfer(const DistrArray &array) {
  Hid plist(array, H5Pcreate(H5P_DATASET_XFER), H5Pclose, "create a transfer property list");
  check(array, H5Pset_dxpl_mpio(plist, H5FD_MPIO_INDEPENDENT), "set independent transfer mode");
  return plist;
}
} // namespace

DistrArrayHDF5::DistrArrayHDF5(const std::shared_ptr<util::PHDF5Handle> &file_handle,
                               std::unique_ptr<Distribution> distribution)
    : DistrArrayDisk(std::move(distribution), file_handle->communicator()), m_file_handle(file_handle) {
  if (m_distribution->border().first != 0)
    DistrArray::error("Distribution of array must start from 0");
  DistrArrayHDF5::open_access();
}

DistrArrayHDF5::DistrArrayHDF5(const std::shared_ptr<util::PHDF5Handle> &file_handle, size_t dimension)
    : DistrArrayHDF5(file_handle, std::make_unique<Distribution>(util::make_distribution_spread_remainder<index_type>(
                                      dimension, mpi_size(file_handle->communicator())))) {}

DistrArrayHDF5::DistrArrayHDF5(const DistrArray &source, const std::shared_ptr<util::PHDF5Handle> &file_handle)
    : DistrArrayHDF5(file_handle, std::make_unique<Distribution>(source.distribution())) {
  DistrArrayHDF5::copy(source);
}

DistrArrayHDF5::DistrArrayHDF5(DistrArrayHDF5 &&source) noexcept
    : DistrArrayDisk(std::move(source)), m_file_handle(std::move(source.m_file_handle)), m_dataset(source.m_dataset) {
  source.m_dataset = DistrArrayHDF5::dataset_default;
}

DistrArrayHDF5 &DistrArrayHDF5::operator=(DistrArrayHDF5 &&source) noexcept {
  DistrArrayHDF5 t{std::move(source)};
  swap(*this, t);
  return *this;
}

DistrArrayHDF5::~DistrArrayHDF5() {
  // Failures are not reported: a destructor must not throw, and the array is gone anyway
  if (dataset_is_open()) {
    H5Dclose(m_dataset);
    m_file_handle->close_file();
  }
}

DistrArrayHDF5 DistrArrayHDF5::CreateTempCopy(const DistrArray &source, const std::string &base_name) {
  const auto name = base_name.empty() ? temp_copy_directory(source, base_name) : base_name;
  return DistrArrayHDF5(source, temp_file_handle(name, source.communicator()));
}

DistrArrayHDF5 DistrArrayHDF5::CreateTempCopy(Span<const value_type> source, MPI_Comm comm,
                                              const std::string &base_name) {
  auto t = DistrArrayHDF5(temp_file_handle(base_name, comm), source.size());
  t.copy(source);
  return t;
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
  swap(x.m_allocated, y.m_allocated);
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
    error("DistrArrayHDF5: failed to open file " + m_file_handle->file_name());
  id = m_file_handle->open_group();
  if (id < 0)
    error("DistrArrayHDF5: failed to open group " + m_file_handle->group_name() + " in file " +
          m_file_handle->file_name());
  auto exists = check(*this, dataset_exists(), "check whether the dataset exists");
  if (exists > 0) {
    if (m_dataset == dataset_default)
      m_dataset =
          check(*this, H5Dopen(m_file_handle->group_id(), dataset_name.c_str(), H5P_DEFAULT), "open the dataset");
  } else {
    hsize_t dimensions[1] = {m_dimension};
    Hid space(*this, H5Screate_simple(1, dimensions, nullptr), H5Sclose, "create a dataspace");
    m_dataset = check(*this,
                      H5Dcreate(m_file_handle->group_id(), dataset_name.c_str(), H5T_NATIVE_DOUBLE, space, H5P_DEFAULT,
                                H5P_DEFAULT, H5P_DEFAULT),
                      "create the dataset");
  }
}

void DistrArrayHDF5::close_access() {
  // Flushing first reports failures to write out cached data, which closing the file would not
  const auto flushed = m_dataset == dataset_default ? 0 : H5Fflush(m_dataset, H5F_SCOPE_LOCAL);
  const auto closed = m_dataset == dataset_default ? 0 : H5Dclose(m_dataset);
  m_dataset = dataset_default;
  m_file_handle->close_file();
  check(*this, flushed, "flush the file");
  check(*this, closed, "close the dataset");
}

void DistrArrayHDF5::erase() {
  open_access();
  {
    Hid lapl_id(*this, H5Dget_access_plist(m_dataset), H5Pclose, "access the property list of the dataset");
    check(*this, H5Ldelete(m_file_handle->group_id(), dataset_name.c_str(), lapl_id), "delete the dataset");
  }
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
  Hid memspace(*this, H5Screate_simple(1, count, nullptr), H5Sclose, "create a dataspace");
  Hid fspace(*this, H5Dget_space(m_dataset), H5Sclose, "get the dataspace of the dataset");
  check(*this, H5Sselect_hyperslab(fspace, H5S_SELECT_SET, offset, nullptr, count, nullptr), "select a hyperslab");
  check(*this, H5Dread(m_dataset, H5T_NATIVE_DOUBLE, memspace, fspace, independent_transfer(*this), buf),
        "read the dataset");
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
  if (lo >= hi)
    return;
  if (!dataset_is_open())
    error("call open_access() before RMA operations");
  hsize_t count[1] = {hi - lo};
  hsize_t offset[1] = {lo};
  Hid memspace(*this, H5Screate_simple(1, count, nullptr), H5Sclose, "create a dataspace");
  Hid fspace(*this, H5Dget_space(m_dataset), H5Sclose, "get the dataspace of the dataset");
  check(*this, H5Sselect_hyperslab(fspace, H5S_SELECT_SET, offset, nullptr, count, nullptr), "select a hyperslab");
  check(*this, H5Dwrite(m_dataset, H5T_NATIVE_DOUBLE, memspace, fspace, independent_transfer(*this), data),
        "write the dataset");
}

void DistrArrayHDF5::acc(index_type lo, index_type hi, const value_type *data) {
  if (lo >= hi)
    return;
  auto disk_copy = get(lo, hi);
  std::transform(begin(disk_copy), end(disk_copy), data, begin(disk_copy), [](auto &l, auto &r) { return l + r; });
  put(lo, hi, &disk_copy[0]);
}

std::vector<DistrArrayHDF5::value_type> DistrArrayHDF5::gather(const std::vector<index_type> &indices) const {
  if (!dataset_is_open())
    error("call open_access() before RMA operations");
  auto sz = indices.size();
  auto data = std::vector<value_type>(sz);
  if (sz == 0)
    return data;
  hsize_t count[] = {sz};
  Hid memspace(*this, H5Screate_simple(1, count, nullptr), H5Sclose, "create a dataspace");
  Hid fspace(*this, H5Dget_space(m_dataset), H5Sclose, "get the dataspace of the dataset");
  auto indices_copy = std::vector<hsize_t>(sz);
  std::copy(begin(indices), end(indices), begin(indices_copy));
  check(*this, H5Sselect_elements(fspace, H5S_SELECT_SET, sz, &indices_copy[0]), "select elements");
  check(*this, H5Dread(m_dataset, H5T_NATIVE_DOUBLE, memspace, fspace, independent_transfer(*this), &data[0]),
        "read the dataset");
  return data;
}

void DistrArrayHDF5::scatter(const std::vector<index_type> &indices, const std::vector<value_type> &data) {
  if (!dataset_is_open())
    error("call open_access() before RMA operations");
  auto sz = indices.size();
  if (sz == 0)
    return;
  hsize_t count[] = {sz};
  Hid memspace(*this, H5Screate_simple(1, count, nullptr), H5Sclose, "create a dataspace");
  Hid fspace(*this, H5Dget_space(m_dataset), H5Sclose, "get the dataspace of the dataset");
  auto indices_copy = std::vector<hsize_t>(sz);
  std::copy(begin(indices), end(indices), begin(indices_copy));
  check(*this, H5Sselect_elements(fspace, H5S_SELECT_SET, sz, &indices_copy[0]), "select elements");
  check(*this, H5Dwrite(m_dataset, H5T_NATIVE_DOUBLE, memspace, fspace, independent_transfer(*this), &data[0]),
        "write the dataset");
}

void DistrArrayHDF5::scatter_acc(std::vector<index_type> &indices, const std::vector<value_type> &data) {
  auto disk_copy = gather(indices);
  std::transform(begin(data), end(data), begin(disk_copy), begin(disk_copy), [](auto &l, auto &r) { return l + r; });
  scatter(indices, disk_copy);
}

std::vector<DistrArrayHDF5::value_type> DistrArrayHDF5::vec() const { return get(0, m_dimension); }

std::shared_ptr<util::PHDF5Handle> DistrArrayHDF5::file_handle() const { return m_file_handle; }

std::string DistrArrayHDF5::file_directory() const {
  if (not m_file_handle)
    return {};
  return fs::absolute(fs::path(m_file_handle->file_name())).parent_path().string();
}

hid_t DistrArrayHDF5::dataset() const { return m_dataset; }

} // namespace molpro::linalg::array
