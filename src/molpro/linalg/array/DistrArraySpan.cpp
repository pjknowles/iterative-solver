#include "DistrArraySpan.h"
#include "util/Distribution.h"
#include <string>
#include <memory>

namespace molpro::linalg::array {

namespace {
int mpi_size(MPI_Comm comm) {
  int rank;
  MPI_Comm_size(comm, &rank);
  return rank;
}
}

DistrArraySpan::DistrArraySpan(size_t dimension, MPI_Comm commun)
    : DistrArraySpan(std::make_unique<Distribution>(
      util::make_distribution_spread_remainder<index_type>(dimension, mpi_size(commun))), commun) {}
      
DistrArraySpan::DistrArraySpan(std::unique_ptr<Distribution> distribution, MPI_Comm commun)
    : DistrArray(distribution->border().second, commun), m_distribution(std::move(distribution)) {
  if (m_distribution->border().first != 0)
    DistrArray::error("Distribution of array must start from 0");
    //TODO: Empty on construction - OK???
}

// DistrArraySpan::~DistrArraySpan() {}

DistrArraySpan::DistrArraySpan(const DistrArraySpan& source)
    : DistrArray(source.size(), source.communicator()),
       m_distribution(source.m_distribution ? std::make_unique<Distribution>(*source.m_distribution) : nullptr) {
  if (!source.empty()) {
    DistrArraySpan::allocate_buffer(source.m_span);
  }
}

DistrArraySpan::DistrArraySpan(const DistrArray& source)
    : DistrArray(source), m_distribution(std::make_unique<Distribution>(source.distribution())) {
  if (!source.empty()) {
    DistrArraySpan::allocate_buffer(Span<value_type>(&(*source.local_buffer())[0], source.size()));
  }
}

DistrArraySpan::DistrArraySpan(DistrArraySpan&& source) noexcept
    : DistrArray(source.m_dimension, source.m_communicator), m_span(std::move(source.m_span)),
      m_allocated(source.m_allocated), m_distribution(std::move(source.m_distribution)) {
  source.m_allocated = false; //TODO: and source.m_span?
}

DistrArraySpan& DistrArraySpan::operator=(const DistrArraySpan& source) {
  if (this == &source)
    return *this;
  if (source.empty() || empty() || !compatible(source)) { // TODO: not sure about this case
    free_buffer();
    DistrArraySpan t{source};
    swap(*this, t);
  } else {
    allocate_buffer(source.m_span);
  }
  return *this;
}

DistrArraySpan& DistrArraySpan::operator=(DistrArraySpan&& source) noexcept {
  DistrArraySpan t{std::move(source)};
  swap(*this, t);
  return *this;
}

void swap(DistrArraySpan& a1, DistrArraySpan& a2) noexcept {
  using std::swap;
  swap(a1.m_dimension, a2.m_dimension);
  swap(a1.m_communicator, a2.m_communicator);
  swap(a1.m_distribution, a2.m_distribution);
  swap(a1.m_allocated, a2.m_allocated);
  swap(a1.m_span, a2.m_span);
}

void DistrArraySpan::free_buffer() {
  if (m_allocated) {
    // TODO: what to do with m_span?
    m_allocated = false;
  }
}

void DistrArraySpan::allocate_buffer() {} //TODO: Doesn't make sense to have it... Or should we allocate empty space?

void DistrArraySpan::allocate_buffer(Span<value_type> buffer) {
  //if (!empty())  // TODO: OK to be "re-writable"?
  //  return;
  if (!m_distribution)
    error("Cannot allocate an array without distribution");
  //m_buffer = std::make_unique<LocalBufferSpan>(buffer);
  m_span = buffer;
  int rank = MPI_Comm_rank(m_communicator, &rank);
  index_type lo, hi;
  std::tie(lo, hi) = m_distribution->range(rank);
  size_t n = hi - lo;
  if (m_span.size() < n)
    error("Specified external buffer is too small");
  m_allocated = true;
}

DistrArraySpan::LocalBufferSpan::LocalBufferSpan(Span<value_type>& source) {
  // TODO: is there a smarter way of doing this (i.e. using Span copy-constructor)???
  m_buffer = source.data();
  m_size = source.size();
}

DistrArraySpan::LocalBufferSpan::LocalBufferSpan(const Span<value_type>& source) {
  // TODO: is there a smarter way of doing this (i.e. using Span copy-constructor)???
  m_buffer = const_cast<value_type*>(source.data());
  m_size = source.size();
}

bool DistrArraySpan::empty() const { return !m_allocated; }

const DistrArray::Distribution& DistrArraySpan::distribution() const {
  if (!m_distribution)
    error("allocate buffer before asking for distribution");
  return *m_distribution;
}

std::unique_ptr<DistrArray::LocalBuffer> DistrArraySpan::local_buffer() {
  return std::make_unique<LocalBufferSpan>(m_span);
}

std::unique_ptr<const DistrArray::LocalBuffer> DistrArraySpan::local_buffer() const {
  return std::make_unique<const LocalBufferSpan>(m_span);
}

DistrArray::value_type DistrArraySpan::at(DistrArray::index_type ind) const {
  value_type val;
  get(ind, ind + 1, &val);
  return val;
}

void DistrArraySpan::set(DistrArray::index_type ind, DistrArray::value_type val) { put(ind, ind + 1, &val); }

void DistrArraySpan::get(DistrArray::index_type lo, DistrArray::index_type hi, DistrArray::value_type* buf) const {
  if (lo >= hi)
    return;
  int rank;
  MPI_Comm_rank(m_communicator, &rank);
  DistrArray::index_type lo_loc, hi_loc;
  std::tie(lo_loc, hi_loc) = m_distribution->range(rank);
  if (lo < lo_loc || hi > hi_loc) {
    error("Only local array indices can be accessed via DistrArraySpan.get() function");
  }
  DistrArray::index_type offset = lo - lo_loc;
  DistrArray::index_type length = hi - lo;
  for (int i = offset; i < offset+length; i++) {
    buf[i-offset] = m_span[i];
  }
}

std::vector<DistrArraySpan::value_type> DistrArraySpan::get(DistrArray::index_type lo,
                                                            DistrArray::index_type hi) const {
  if (lo >= hi)
    return {};
  auto buf = std::vector<DistrArray::value_type>(hi - lo);
  get(lo, hi, &buf[0]);
  return buf;
}

void DistrArraySpan::put(DistrArray::index_type lo, DistrArray::index_type hi, const DistrArray::value_type* data) {
  if (lo >= hi)
    return;
  int rank;
  MPI_Comm_rank(m_communicator, &rank);
  DistrArray::index_type lo_loc, hi_loc;
  std::tie(lo_loc, hi_loc) = m_distribution->range(rank);
  if (lo < lo_loc || hi > hi_loc) {
    error("Only values at local array indices can be written via DistrArraySpan.put() function");
  }
  DistrArray::index_type offset = lo - lo_loc;
  DistrArray::index_type length = hi - lo;
  for (int i = offset; i < offset+length; i++) {
    m_span[i] = data[i-offset];
  }
}

void DistrArraySpan::acc(DistrArray::index_type lo, DistrArray::index_type hi, const DistrArray::value_type* data) {
  if (lo >= hi)
    return;
  auto disk_copy = get(lo, hi);
  std::transform(disk_copy.begin(), disk_copy.end(), data, disk_copy.begin(), [](auto& l, auto& r) { return l + r; });
  put(lo, hi, &disk_copy[0]);
}

std::vector<DistrArraySpan::value_type> DistrArraySpan::gather(const std::vector<index_type>& indices) const {
  std::vector<value_type> data;
  data.reserve(indices.size());
  int rank;
  MPI_Comm_rank(m_communicator, &rank);
  auto minmax = std::minmax_element(indices.begin(), indices.end());
  DistrArray::index_type lo_loc, hi_loc;
  std::tie(lo_loc, hi_loc) = m_distribution->range(rank);
  if (*minmax.first < lo_loc || *minmax.second > hi_loc) {
    error("Only local array indices can be accessed via DistrArraySpan.gather() function");
  }
  for (auto i : indices) {
    data.push_back(at(i));
  }
  return data;
}

void DistrArraySpan::scatter(const std::vector<index_type>& indices, const std::vector<value_type>& data) {
  if (indices.size() != data.size()) {
    error("Length of the indices and data vectors should be the same: DistrArray::scatter()");
  }
  int rank;
  MPI_Comm_rank(m_communicator, &rank);
  auto minmax = std::minmax_element(indices.begin(), indices.end());
  DistrArray::index_type lo_loc, hi_loc;
  std::tie(lo_loc, hi_loc) = m_distribution->range(rank);
  if (*minmax.first < lo_loc || *minmax.second > hi_loc) {
    error("Only local array indices can be accessed via DistrArraySpan.gather() function");
  }
  for (auto i : indices) {
    set(i, data[i-*minmax.first]);
  }
}

void DistrArraySpan::scatter_acc(std::vector<index_type>& indices, const std::vector<value_type>& data) {
  auto disk_copy = gather(indices);
  std::transform(data.begin(), data.end(), disk_copy.begin(), disk_copy.begin(),
                 [](auto& l, auto& r) { return l + r; });
  scatter(indices, disk_copy);
}

std::vector<DistrArraySpan::value_type> DistrArraySpan::vec() const {
  int rank;
  MPI_Comm_rank(m_communicator, &rank);
  DistrArray::index_type lo_loc, hi_loc;
  std::tie(lo_loc, hi_loc) = m_distribution->range(rank);
  return get(lo_loc, hi_loc);
}

} // molpro::linalg::array