#include "DistrArrayMPI3.h"
#include "util.h"
#include "util/Distribution.h"
#include <algorithm>
#include <tuple>
namespace molpro {
namespace linalg {
namespace array {

namespace {
int comm_size(MPI_Comm comm) {
  int res;
  MPI_Comm_size(comm, &res);
  return res;
}

int comm_rank(MPI_Comm comm) {
  int res;
  MPI_Comm_rank(comm, &res);
  return res;
}
} // namespace

DistrArrayMPI3::DistrArrayMPI3() = default;

DistrArrayMPI3::DistrArrayMPI3(size_t dimension, MPI_Comm commun, std::shared_ptr<Profiler> prof)
    : DistrArrayMPI3(std::make_unique<Distribution>(
                         util::make_distribution_spread_remainder<index_type>(dimension, comm_size(commun))),
                     commun, std::move(prof)) {}

DistrArrayMPI3::DistrArrayMPI3(std::unique_ptr<Distribution> distribution, MPI_Comm commun,
                               std::shared_ptr<Profiler> prof)
    : DistrArray(distribution->border().second, commun, std::move(prof)), m_distribution(std::move(distribution)) {
  if (m_distribution->border().first != 0)
    DistrArray::error("Distribution of array must start from 0");
}

DistrArrayMPI3::~DistrArrayMPI3() {
  if (m_allocated)
    DistrArrayMPI3::free_buffer();
}

void DistrArrayMPI3::allocate_buffer() {
  if (!empty())
    return;
  int rank;
  MPI_Comm_rank(m_communicator, &rank);
  index_type lo, hi;
  std::tie(lo, hi) = m_distribution->range(rank);
  MPI_Aint n = hi - lo;
  double* base = nullptr;
  int size_of_type = sizeof(value_type);
  n *= size_of_type;
  MPI_Win_allocate(n, size_of_type, MPI_INFO_NULL, m_communicator, &base, &m_win);
  MPI_Win_lock_all(0, m_win);
  m_allocated = true;
}

bool DistrArrayMPI3::empty() const { return !m_allocated; }

DistrArrayMPI3::DistrArrayMPI3(const DistrArrayMPI3& source)
    : DistrArray(source.size(), source.communicator(), source.m_prof),
      m_distribution(source.m_distribution ? std::make_unique<Distribution>(*source.m_distribution) : nullptr) {
  if (!source.empty()) {
    DistrArrayMPI3::allocate_buffer();
    DistrArray::copy(source);
  }
}

DistrArrayMPI3::DistrArrayMPI3(DistrArrayMPI3&& source) noexcept
    : DistrArray(source.m_dimension, source.m_communicator, source.m_prof), m_win(source.m_win),
      m_allocated(source.m_allocated), m_distribution(std::move(source.m_distribution)) {
  source.m_allocated = false;
}

DistrArrayMPI3& DistrArrayMPI3::operator=(const DistrArrayMPI3& source) {
  if (this == &source)
    return *this;
  if (source.empty() || empty() || !compatible(source)) {
    free_buffer();
    DistrArrayMPI3 t{source};
    swap(*this, t);
  } else {
    allocate_buffer();
    copy(source);
    m_prof = source.m_prof;
  }
  return *this;
}

DistrArrayMPI3& DistrArrayMPI3::operator=(DistrArrayMPI3&& source) noexcept {
  DistrArrayMPI3 t{std::move(source)};
  swap(*this, t);
  return *this;
}

void swap(DistrArrayMPI3& a1, DistrArrayMPI3& a2) noexcept {
  using std::swap;
  swap(a1.m_dimension, a2.m_dimension);
  swap(a1.m_communicator, a2.m_communicator);
  swap(a1.m_prof, a2.m_prof);
  swap(a1.m_distribution, a2.m_distribution);
  swap(a1.m_allocated, a2.m_allocated);
  swap(a1.m_win, a2.m_win);
}

void DistrArrayMPI3::free_buffer() {
  if (m_allocated) {
    MPI_Win_unlock_all(m_win);
    MPI_Win_free(&m_win);
    m_allocated = false;
  }
}

void DistrArrayMPI3::sync() const {
  if (!empty()) {
    MPI_Win_flush_all(m_win);
    MPI_Win_sync(m_win);
  }
  MPI_Barrier(m_communicator);
}

std::unique_ptr<const DistrArray::LocalBuffer> DistrArrayMPI3::local_buffer() const {
  return std::make_unique<const LocalBufferMPI3>(*const_cast<DistrArrayMPI3*>(this));
}

std::unique_ptr<DistrArrayMPI3::LocalBuffer> DistrArrayMPI3::local_buffer() {
  return std::make_unique<LocalBufferMPI3>(*this);
}

DistrArray::value_type DistrArrayMPI3::at(index_type ind) const {
  value_type val;
  get(ind, ind, &val);
  return val;
}
void DistrArrayMPI3::set(index_type ind, value_type val) { put(ind, ind, &val); }

void DistrArrayMPI3::_get_put(index_type lo, index_type hi, const value_type* buf, RMAType option) {
  if (lo > hi)
    return;
  auto name = std::string{"DistrArrayMPI3::_get_put"};
  if (hi >= m_dimension)
    error(name + " out of bounds");
  if (empty())
    error(name + " called on an empty array");
  util::ScopeProfiler prof{m_prof, name};
  index_type p_lo, p_hi;
  std::tie(p_lo, p_hi) = m_distribution->cover(lo, hi);
  auto* curr_buf = const_cast<value_type*>(buf);
  auto requests = std::vector<MPI_Request>(p_hi - p_lo + 1);
  for (size_t i = p_lo; i < p_hi + 1; ++i) {
    index_type bound_lo, bound_hi;
    std::tie(bound_lo, bound_hi) = m_distribution->range(i);
    --bound_hi;
    auto local_lo = std::max(lo, bound_lo);
    auto local_hi = std::min(hi, bound_hi);
    MPI_Aint offset = (local_lo - bound_lo);
    int count = (int(local_hi - local_lo) + 1);
    if (option == RMAType::get)
      MPI_Rget(curr_buf, count, MPI_DOUBLE, i, offset, count, MPI_DOUBLE, m_win, &requests[i - p_lo]);
    else if (option == RMAType::put)
      MPI_Rput(curr_buf, count, MPI_DOUBLE, i, offset, count, MPI_DOUBLE, m_win, &requests[i - p_lo]);
    else if (option == RMAType::acc)
      MPI_Raccumulate(curr_buf, count, MPI_DOUBLE, i, offset, count, MPI_DOUBLE, MPI_SUM, m_win, &requests[i - p_lo]);
    curr_buf += count;
  }
  MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
}

void DistrArrayMPI3::get(index_type lo, index_type hi, value_type* buf) const {
  const_cast<DistrArrayMPI3*>(this)->_get_put(lo, hi, buf, RMAType::get);
}

std::vector<DistrArrayMPI3::value_type> DistrArrayMPI3::get(index_type lo, index_type hi) const {
  if (lo > hi)
    return {};
  auto val = std::vector<value_type>(hi - lo + 1);
  get(lo, hi, val.data());
  return val;
}

void DistrArrayMPI3::put(index_type lo, index_type hi, const value_type* data) { _get_put(lo, hi, data, RMAType::put); }

std::vector<DistrArrayMPI3::value_type> DistrArrayMPI3::gather(const std::vector<index_type>& indices) const {
  auto data = std::vector<value_type>(indices.size());
  const_cast<DistrArrayMPI3*>(this)->_gather_scatter(indices, data, RMAType::gather);
  return data;
}

void DistrArrayMPI3::scatter(const std::vector<index_type>& indices, const std::vector<value_type>& data) {
  _gather_scatter(indices, const_cast<std::vector<value_type>&>(data), RMAType::scatter);
}

void DistrArrayMPI3::scatter_acc(std::vector<index_type>& indices, const std::vector<value_type>& data) {
  _gather_scatter(indices, const_cast<std::vector<value_type>&>(data), RMAType::scatter_acc);
}

void DistrArrayMPI3::_gather_scatter(const std::vector<index_type>& indices, std::vector<value_type>& data,
                                     RMAType option) {
  if (indices.empty())
    return;
  auto name = std::string{"DistrArrayMPI3::_gather_scatter"};
  if (*std::max_element(indices.begin(), indices.end()) >= m_dimension)
    error(name + " out of bounds");
  if (indices.size() > data.size())
    error(name + " data buffer is too small");
  if (empty())
    error(name + " called on an empty array");
  util::ScopeProfiler prof{m_prof, name};
  auto requests = std::vector<MPI_Request>(indices.size());
  for (size_t i = 0; i < indices.size(); ++i) {
    int p;
    index_type lo;
    p = m_distribution->cover(indices[i]);
    std::tie(lo, std::ignore) = m_distribution->range(p);
    MPI_Aint offset = indices[i] - lo;
    if (option == RMAType::gather)
      MPI_Rget(&data[i], 1, MPI_DOUBLE, p, offset, 1, MPI_DOUBLE, m_win, &requests[i]);
    else if (option == RMAType::scatter)
      MPI_Rput(&data[i], 1, MPI_DOUBLE, p, offset, 1, MPI_DOUBLE, m_win, &requests[i]);
    else if (option == RMAType::scatter_acc)
      MPI_Raccumulate(&data[i], 1, MPI_DOUBLE, p, offset, 1, MPI_DOUBLE, MPI_SUM, m_win, &requests[i]);
  }
  MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
}
std::vector<DistrArrayMPI3::value_type> DistrArrayMPI3::vec() const { return get(0, m_dimension - 1); }

void DistrArrayMPI3::acc(index_type lo, index_type hi, const value_type* data) { _get_put(lo, hi, data, RMAType::acc); }

void DistrArrayMPI3::error(const std::string& message) const { MPI_Abort(m_communicator, 1); }

const DistrArray::Distribution& DistrArrayMPI3::distribution() const {
  if (!m_distribution)
    error("allocate buffer before asking for distribution");
  return *m_distribution;
}

DistrArrayMPI3::LocalBufferMPI3::LocalBufferMPI3(DistrArrayMPI3& source) {
  if (!source.m_allocated)
    source.error("attempting to access local buffer of empty array");
  int rank;
  MPI_Comm_rank(source.communicator(), &rank);
  index_type hi;
  std::tie(m_start, hi) = source.distribution().range(rank);
  m_size = hi - m_start;
  int flag;
  MPI_Win_get_attr(source.m_win, MPI_WIN_BASE, &m_buffer, &flag);
}

} // namespace array
} // namespace linalg
} // namespace molpro
