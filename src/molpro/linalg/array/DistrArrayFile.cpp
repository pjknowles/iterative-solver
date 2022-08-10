#include <unistd.h>

#include <cstring>
#include <molpro/linalg/array/DistrArrayFile.h>
#include <molpro/linalg/array/util/Distribution.h>
#include <molpro/linalg/array/util/temp_file.h>
#include <utility>

#include <memory>
#include <molpro/Profiler.h>
#ifdef _WIN32
#include <stdio.h>
#else
#include <sys/resource.h>
#include <sys/time.h>
#include <sys/types.h>
#endif
#include <iostream>
#include <mutex>

namespace molpro::linalg::array {
namespace {
int mpi_size(MPI_Comm comm) {
  if (comm == molpro::mpi::comm_global())
    return molpro::mpi::size_global();
#ifdef HAVE_MPI_H
  int size;
  MPI_Comm_size(comm, &size);
  return size;
#endif
  throw std::logic_error("Attempt to access MPI communicator in serial mode");
}

int mpi_rank(MPI_Comm comm) {
  if (comm == molpro::mpi::comm_global())
    return molpro::mpi::rank_global();
#ifdef HAVE_MPI_H
  int rank;
  MPI_Comm_rank(comm, &rank);
  return rank;
#endif
  throw std::logic_error("Attempt to access MPI communicator in serial mode");
}
} // namespace

std::tuple<DistrArrayFile::index_type, DistrArrayFile::index_type, DistrArrayFile::index_type>
DistrArrayFile::local_bounds() const {
  index_type lo_loc, hi_loc, size_loc;
  std::tie(lo_loc, hi_loc) = m_distribution->range(mpi_rank(m_communicator));
  size_loc = hi_loc - lo_loc;
  return {lo_loc, hi_loc, size_loc};
}

DistrArrayFile::DistrArrayFile(size_t dimension, MPI_Comm comm, const std::string& directory)
    : DistrArrayFile(std::make_unique<Distribution>(
                         util::make_distribution_spread_remainder<index_type>(dimension, mpi_size(comm))),
                     comm, directory) {}

std::mutex s_open_error_mutex;
DistrArrayFile::DistrArrayFile(std::unique_ptr<Distribution> distribution, MPI_Comm comm, const std::string& directory)
    : DistrArrayDisk(std::move(distribution), comm), m_directory(fs::absolute(fs::path(directory))),
      m_filename(util::temp_file_name((fs::path(m_directory) / "DistrArrayFile-"), ".dat")),
      m_stream(std::make_unique<std::fstream>(m_filename.c_str(),
                                              std::ios::out | std::ios::binary | std::ios::trunc | std::ios::in)) {
  {
    auto prof = molpro::Profiler::single()->push("DistrArrayFile()");
    prof++;
#ifdef _WIN32
    _setmaxstdio(8192);
#else
    struct rlimit rlim;
    getrlimit(RLIMIT_NOFILE, &rlim);
    if (rlim.rlim_cur < rlim.rlim_max) {
      rlim.rlim_cur = rlim.rlim_max;
      setrlimit(RLIMIT_NOFILE, &rlim);
    }
#endif
  }
  if (m_distribution->border().first != 0)
    DistrArray::error("Distribution of array must start from 0");
  if (not m_stream->is_open()) {
    std::lock_guard<std::mutex> lock(s_open_error_mutex);
    std::cerr << "Failure to open file " << m_filename << ": " << std::strerror(errno) << std::endl;
    DistrArray::error(std::string{"Failure to open file, "} + std::strerror(errno));
  }
#if !defined(_WIN32) && !defined(WIN32)
  fs::remove(m_filename);
#endif
}

DistrArrayFile::DistrArrayFile(const DistrArrayFile& source)
    : DistrArrayFile(std::make_unique<Distribution>(source.distribution()), source.communicator()) {
  DistrArrayFile::copy(source);
}

DistrArrayFile::DistrArrayFile(DistrArrayFile&& source)
    : DistrArrayDisk(std::move(source)), m_directory(std::move(source.m_directory)),
      m_filename(std::move(source.m_filename)), m_stream(std::move(source.m_stream)) {}

DistrArrayFile::DistrArrayFile(const DistrArray& source)
    : DistrArrayFile(std::make_unique<Distribution>(source.distribution()), source.communicator()) {
  DistrArrayFile::copy(source);
}

DistrArrayFile& DistrArrayFile::operator=(DistrArrayFile&& source) noexcept {
  DistrArrayFile t{std::move(source)};
  swap(*this, t);
  return *this;
}

DistrArrayFile DistrArrayFile::CreateTempCopy(const DistrArray& source, const std::string& directory) {
  DistrArrayFile t(std::make_unique<Distribution>(source.distribution()), source.communicator(), directory);
  t.copy(source);
  return t;
}

void swap(DistrArrayFile& x, DistrArrayFile& y) noexcept {
  using std::swap;
  swap(x.m_dimension, y.m_dimension);
  swap(x.m_communicator, y.m_communicator);
  swap(x.m_allocated, y.m_allocated);
  swap(x.m_distribution, y.m_distribution);
  swap(x.m_directory, y.m_directory);
  swap(x.m_filename, y.m_filename);
  swap(x.m_stream, y.m_stream);
}

DistrArrayFile::~DistrArrayFile() {
  auto prof = molpro::Profiler::single()->push("~DistrArrayFile()");
  prof++;
  if (m_stream.get() != nullptr) {
    m_stream->close();
//    std::cout << "closing " << m_filename << std::endl;
    if (fs::exists(m_filename))
      fs::remove(m_filename);
    m_stream.release();
  }
}

bool DistrArrayFile::compatible(const DistrArrayFile& source) const {
  auto res = DistrArray::compatible(source);
  if (m_distribution && source.m_distribution)
    res &= m_distribution->compatible(*source.m_distribution);
  else
    res &= !m_distribution && !source.m_distribution;
  return res;
}

DistrArray::value_type DistrArrayFile::at(DistrArray::index_type ind) const {
  value_type val;
  get(ind, ind + 1, &val);
  return val;
}

void DistrArrayFile::set(DistrArray::index_type ind, DistrArray::value_type val) { put(ind, ind + 1, &val); }

void DistrArrayFile::get(DistrArray::index_type lo, DistrArray::index_type hi, DistrArray::value_type* buf) const {
  auto lock = std::lock_guard<std::mutex>(m_mutex);
  if (lo >= hi)
    return;
  DistrArray::index_type length = hi - lo;
  DistrArray::index_type lo_loc, hi_loc;
  auto bounds_loc = local_bounds();
  std::tie(lo_loc, hi_loc) = {std::get<0>(bounds_loc), std::get<1>(bounds_loc)};
  if (lo < lo_loc || hi > hi_loc) {
    error(
        "Only local array indices can be accessed via DistrArrayFile.get() function\nlo:lo_loc=" + std::to_string(lo) +
        ":" + std::to_string(lo_loc) + "\nhi:hi_loc=" + std::to_string(hi) + ":" + std::to_string(hi_loc));
  }
  DistrArray::index_type offset = lo - lo_loc;
  m_stream->seekg(0, std::ios::end);
  const auto wordlength = sizeof(DistrArray::value_type);
  size_t current = m_stream->tellg() / wordlength;
  if (current < (offset + length)) {
    std::fill(buf + current - offset, buf + length, 0);
  }
  m_stream->seekg(offset * wordlength);
  const auto readlength = current - offset > 0 ? std::min(length, current - offset) * wordlength : 0;
//  if (readlength==0)
//  std::cout << "current=" << current << ", offset=" << offset << ", length=" << length << ", readlength=" << readlength
//            << std::endl;
  if (readlength == 0)
    return;
  m_stream->read((char*)buf, readlength);
  if (m_stream->fail()) {
    std::lock_guard<std::mutex> lock(s_open_error_mutex);
    std::cerr << "Error in reading file " << m_filename << ": " << std::strerror(errno) << std::endl;
    throw std::runtime_error(std::string{"Error in reading ,"} + std::strerror(errno));
  }
}

std::vector<DistrArrayFile::value_type> DistrArrayFile::get(DistrArray::index_type lo,
                                                            DistrArray::index_type hi) const {
  //  std::cout << "DistrArrayFile::get() " << m_filename<<" "<<lo<<" "<<hi<<std::endl;
  if (lo >= hi)
    return {};
  auto buf = std::vector<DistrArray::value_type>(hi - lo);
  get(lo, hi, &buf[0]);
  return buf;
}

void DistrArrayFile::put(DistrArray::index_type lo, DistrArray::index_type hi, const DistrArray::value_type* data) {
  //  std::cout << "DistrArrayFile::put() " << m_filename<<" "<<lo<<" "<<hi<<std::endl;
  auto lock = std::lock_guard<std::mutex>(m_mutex);
  if (lo >= hi)
    return;
  auto prof = molpro::Profiler::single()->push("DistrArrayFile::put()");
  prof += hi - lo;
  DistrArray::index_type lo_loc, hi_loc;
  auto bounds_loc = local_bounds();
  std::tie(lo_loc, hi_loc) = {std::get<0>(bounds_loc), std::get<1>(bounds_loc)};
  if (lo < lo_loc || hi > hi_loc) {
    error("Only values at local array indices can be written via DistrArrayFile.put() function");
  }
  DistrArray::index_type offset = lo - lo_loc;
  DistrArray::index_type length = hi - lo;
  m_stream->seekp(offset * sizeof(DistrArray::value_type));
  m_stream->write((const char*)data, length * sizeof(DistrArray::value_type));
}

void DistrArrayFile::acc(DistrArray::index_type lo, DistrArray::index_type hi, const DistrArray::value_type* data) {
  if (lo >= hi)
    return;
  auto prof = molpro::Profiler::single()->push("DistrArrayFile::acc()");
  auto disk_copy = get(lo, hi);
  prof += hi - lo;
  std::transform(disk_copy.begin(), disk_copy.end(), data, disk_copy.begin(), [](auto& l, auto& r) { return l + r; });
  put(lo, hi, &disk_copy[0]);
}

std::vector<DistrArrayFile::value_type> DistrArrayFile::gather(const std::vector<index_type>& indices) const {
  std::vector<value_type> data;
  data.reserve(indices.size());
  auto minmax = std::minmax_element(indices.begin(), indices.end());
  DistrArray::index_type lo_loc, hi_loc;
  auto bounds_loc = local_bounds();
  std::tie(lo_loc, hi_loc) = {std::get<0>(bounds_loc), std::get<1>(bounds_loc)};
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
  auto prof = molpro::Profiler::single()->push("DistrArrayFile::scatter()");
  auto minmax = std::minmax_element(indices.begin(), indices.end());
  DistrArray::index_type lo_loc, hi_loc;
  auto bounds_loc = local_bounds();
  std::tie(lo_loc, hi_loc) = {std::get<0>(bounds_loc), std::get<1>(bounds_loc)};
  if (*minmax.first < lo_loc || *minmax.second > hi_loc) {
    error("Only local array indices can be accessed via DistrArrayFile.gather() function");
  }
  for (auto i : indices) {
    set(i, data[i - *minmax.first]); // TODO: check it shouldn't be data[*minmax.first++]??
  }
  prof += indices.size();
}

void DistrArrayFile::scatter_acc(std::vector<index_type>& indices, const std::vector<value_type>& data) {
  auto disk_copy = gather(indices);
  std::transform(data.begin(), data.end(), disk_copy.begin(), disk_copy.begin(),
                 [](auto& l, auto& r) { return l + r; });
  scatter(indices, disk_copy);
}

std::vector<DistrArrayFile::value_type> DistrArrayFile::vec() const {
  DistrArray::index_type lo_loc, hi_loc;
  auto bounds_loc = local_bounds();
  std::tie(lo_loc, hi_loc) = {std::get<0>(bounds_loc), std::get<1>(bounds_loc)};
  return get(lo_loc, hi_loc);
}

} // namespace molpro::linalg::array
