#include <unistd.h>

#include <utility>
#include "DistrArrayFile.h"
#include "util/Distribution.h"
#include "util/temp_file.h"

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

DistrArrayFile::DistrArrayFile() = default;

DistrArrayFile::DistrArrayFile(DistrArrayFile&& source) noexcept : DistrArrayDisk(std::move(source)),
                m_file_name(source.m_file_name), m_file(std::move(source.m_file)) {}
                
DistrArrayFile::DistrArrayFile(size_t dimension, MPI_Comm comm, const std::string &base_name) :
                DistrArrayFile(util::temp_file_name(base_name, ""), dimension, comm){}
                
DistrArrayFile::DistrArrayFile(std::unique_ptr<Distribution> distribution, MPI_Comm comm,
                const std::string& base_name) :
                DistrArrayFile(util::temp_file_name(base_name, ""), std::move(distribution), comm) {}
                
DistrArrayFile::DistrArrayFile(std::string file_name,size_t dimension, MPI_Comm comm) :
                DistrArrayFile(std::move(file_name),
                std::make_unique<Distribution>(util::make_distribution_spread_remainder<index_type>(
                dimension, mpi_size(comm))),comm) {}
                    
DistrArrayFile::DistrArrayFile(std::string file_name, std::unique_ptr<Distribution> distribution, MPI_Comm comm) :
                DistrArrayDisk(std::move(distribution), comm), m_file_name(std::move(file_name)), m_file(make_file()) {
  if (m_distribution->border().first != 0)
    DistrArray::error("Distribution of array must start from 0");
}

DistrArrayFile::DistrArrayFile(std::string file_name, MPI_Comm comm) :
                DistrArrayDisk(std::make_unique<Distribution>(std::vector<index_type>{0, 0}), comm),
                m_file_name(std::move(file_name)), m_file(make_file()) {}

DistrArrayFile::DistrArrayFile(const DistrArray& source) :
                DistrArrayDisk(std::make_unique<Distribution>(source.distribution()), source.communicator()) {
  if (!source.empty()){
    DistrArrayFile::open_access();
    DistrArrayFile::copy(source);
  } else {
    m_file = make_file();
  }
}

DistrArrayFile& DistrArrayFile::operator=(DistrArrayFile&& source) noexcept{
  DistrArrayFile t{std::move(source)};
  swap(*this, t);
  return *this;
}

void swap(DistrArrayFile& x, DistrArrayFile& y) noexcept{
  using std::swap;
  swap(x.m_dimension, y.m_dimension);
  swap(x.m_communicator, y.m_communicator);
  swap(x.m_allocated, y.m_allocated);
  swap(x.m_view_buffer, y.m_view_buffer);
  swap(x.m_owned_buffer, y.m_owned_buffer);
  swap(x.m_distribution, y.m_distribution);
  swap(x.m_file_name, y.m_file_name);
  swap(x.m_file, y.m_file);
}

DistrArrayFile::~DistrArrayFile(){
  // TODO: any point flushing buffer?
  if (m_file.is_open()) {
    DistrArrayFile::close_access();
  }
}

bool DistrArrayFile::compatible(const DistrArrayFile& source) const{
  auto res = DistrArray::compatible(source);
  if (m_distribution && source.m_distribution)
    res &= m_distribution->compatible(*source.m_distribution);
  else
    res &= !m_distribution && !source.m_distribution;
  return res;
}

std::fstream DistrArrayFile::make_file() {
  std::fstream file;
  if (m_file_name.empty()) {
    m_file_name = util::temp_file_name(".temp_array", "");
  }
  file.open(m_file_name.c_str(), std::ios::out | std::ios::binary);
  file.close();
  file.open(m_file_name.c_str(), std::ios::out | std::ios::in | std::ios::binary);
  unlink(m_file_name.c_str());
  return file;
}

void DistrArrayFile::open_access() {
  if (!m_file.is_open()) {
    m_file = make_file();
  }
}
void DistrArrayFile::close_access() {
  if (!m_file.is_open()) {
    error("must provide an existing file stream object before closing access to file disk array");
  } else {
    m_file.close();
  }
}

bool DistrArrayFile::empty() const {
  auto current = m_file.tellg(); // get current position
  m_file.seekg(0, std::ios::end);
  auto res = DistrArrayDisk::empty() && m_file.tellg() == 0;
  m_file.seekg(current, std::ios::beg); // return back to original position
  return res;
}

void DistrArrayFile::erase() {
  // Remove file and then re-create it
  close_access();
  open_access();
}

DistrArray::value_type DistrArrayFile::at(DistrArray::index_type ind) const {
  value_type val;
  get(ind, ind + 1, &val);
  return val;
}

void DistrArrayFile::set(DistrArray::index_type ind, DistrArray::value_type val) {
  put(ind, ind+1, &val);
}

void DistrArrayFile::get(DistrArray::index_type lo, DistrArray::index_type hi, DistrArray::value_type* buf) const {
  if (lo >= hi)
    return;
  int rank;
  MPI_Comm_rank(m_communicator, &rank);
  DistrArray::index_type lo_loc, hi_loc;
  std::tie(lo_loc, hi_loc) = m_distribution->range(rank);
  if (lo < lo_loc || hi > hi_loc) {
    error("Only local array indices can be accessed via DistrArrayFile.get() function");
  }
  DistrArray::index_type offset = lo - lo_loc;
  DistrArray::index_type length = hi - lo;
  m_file.seekg(offset*sizeof(DistrArray::value_type));
  m_file.read((char*)buf, length*sizeof(DistrArray::value_type));
}

std::vector<DistrArrayFile::value_type> DistrArrayFile::get(DistrArray::index_type lo, DistrArray::index_type hi) const {
  if (lo >= hi)
    return {};
  auto buf = std::vector<DistrArray::value_type>(hi - lo);
  get(lo, hi, &buf[0]);
  return buf;
}

void DistrArrayFile::put(DistrArray::index_type lo, DistrArray::index_type hi, const DistrArray::value_type* data) {
  if (lo >= hi)
    return;
  int rank;
  MPI_Comm_rank(m_communicator, &rank);
  DistrArray::index_type lo_loc, hi_loc;
  std::tie(lo_loc, hi_loc) = m_distribution->range(rank);
  if (lo < lo_loc || hi > hi_loc) {
    error("Only values at local array indices can be written via DistrArrayFile.put() function");
  }
  DistrArray::index_type offset = lo - lo_loc;
  DistrArray::index_type length = hi - lo;
  m_file.seekp(offset*sizeof(DistrArray::value_type));
  m_file.write((const char*)data, length*sizeof(DistrArray::value_type));
}

void DistrArrayFile::acc(DistrArray::index_type lo, DistrArray::index_type hi, const DistrArray::value_type* data) {
  if (lo >= hi)
    return;
  auto disk_copy = get(lo, hi);
  std::transform(begin(disk_copy), end(disk_copy), data, begin(disk_copy), [](auto &l, auto &r) { return l + r; });
  put(lo, hi, &disk_copy[0]);
}

std::vector<DistrArrayFile::value_type> DistrArrayFile::gather(const std::vector<index_type>& indices) const {
  auto sz = indices.size();
  auto data = std::vector<value_type>(sz);
  int rank;
  MPI_Comm_rank(m_communicator, &rank);
  auto minmax = std::minmax_element(indices.begin(), indices.end());
  DistrArray::index_type lo_loc, hi_loc;
  std::tie(lo_loc, hi_loc) = m_distribution->range(rank);
  if (*minmax.first < lo_loc || *minmax.second > hi_loc) {
    error("Only local array indices can be accessed via DistrArrayFile.gather() function");
  }
  for (auto i : indices) {
    data.push_back(at(i));
  }
  return data;
}

void DistrArrayFile::scatter(const std::vector<index_type>& indices, const std::vector<value_type>& data) {
  if (indices.size() != data.size()) {
    error("Length of the indices and data vectors should be the same: DistrArray::scatter()");
  }
  int rank;
  MPI_Comm_rank(m_communicator, &rank);
  auto minmax = std::minmax_element(indices.begin(), indices.end());
  DistrArray::index_type lo_loc, hi_loc;
  std::tie(lo_loc, hi_loc) = m_distribution->range(rank);
  if (*minmax.first < lo_loc || *minmax.second > hi_loc) {
    error("Only local array indices can be accessed via DistrArrayFile.gather() function");
  }
  for (auto i : indices) {
    set(i, data[i]);
  }
}

void DistrArrayFile::scatter_acc(std::vector<index_type>& indices, const std::vector<value_type>& data) {
  auto disk_copy = gather(indices);
  std::transform(begin(data), end(data), begin(disk_copy), begin(disk_copy), [](auto &l, auto &r) { return l + r; });
  scatter(indices, disk_copy);
}

std::vector<DistrArrayFile::value_type> DistrArrayFile::vec() const { return get(0, m_dimension); }

std::string DistrArrayFile::file_name() const { return m_file_name; }

} // namespace array
} // namespace linalg
} // namespace molpro
