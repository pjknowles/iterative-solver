#include "DistrArray.h"
#include "util.h"
#include <functional>
#include <molpro/Profiler.h>
#include <numeric>

namespace molpro::gci::array {

DistrArray::DistrArray(size_t dimension, MPI_Comm commun, std::shared_ptr<molpro::Profiler> prof)
    : m_dimension(dimension), m_communicator(commun), m_prof(std::move(prof)) {}

void DistrArray::sync() const { MPI_Barrier(m_communicator); }

void DistrArray::error(const std::string& message) const {
  std::cout << message << std::endl;
  MPI_Abort(m_communicator, 1);
}

size_t DistrArray::size() const { return m_dimension; }

bool DistrArray::compatible(const DistrArray& other) const {
  int comp;
  MPI_Comm_compare(m_communicator, other.m_communicator, &comp);
  return (m_dimension == other.m_dimension) && (MPI_IDENT || MPI_CONGRUENT);
}

bool DistrArray::empty() const { return true; }

void DistrArray::zero() { set(0); }

void DistrArray::set(DistrArray::value_type val) {
  auto p = util::ScopeProfiler(m_prof, "Array::set");
  for (auto& el : *local_buffer())
    el = val;
}

void DistrArray::axpy(DistrArray::value_type a, const DistrArray& x) {
  auto name = std::string{"Array::axpy"};
  if (!compatible(x))
    error(name + " incompatible arrays");
  if (a == 0)
    return;
  auto p = util::ScopeProfiler(m_prof, name);
  auto loc = local_buffer();
  auto loc_x = x.local_buffer();
  if (!loc->compatible(*loc_x))
    error(name + " incompatible local buffers");
  if (a == 1)
    for (size_t i = 0; i < loc->size(); ++i)
      loc->at(i) += loc_x->at(i);
  else if (a == -1)
    for (size_t i = 0; i < loc->size(); ++i)
      loc->at(i) -= loc_x->at(i);
  else
    for (size_t i = 0; i < loc->size(); ++i)
      loc->at(i) += a * loc_x->at(i);
}

void DistrArray::scal(DistrArray::value_type a) {
  auto p = util::ScopeProfiler(m_prof, "Array::scal");
  for (auto& el : *local_buffer())
    el *= a;
}
void DistrArray::add(const DistrArray& x) { axpy(1, x); }

void DistrArray::add(DistrArray::value_type a) {
  auto p = util::ScopeProfiler(m_prof, "Array::add");
  for (auto& el : *local_buffer())
    el += a;
}

void DistrArray::sub(const DistrArray& x) { axpy(-1, x); }

void DistrArray::sub(DistrArray::value_type a) { add(-a); }

void DistrArray::recip() {
  auto p = util::ScopeProfiler(m_prof, "Array::recip");
  for (auto& el : *local_buffer())
    el = 1. / el;
}

void DistrArray::times(const DistrArray& x) {
  auto name = std::string{"Array::times"};
  if (!compatible(x))
    error(name + " incompatible arrays");
  auto p = util::ScopeProfiler(m_prof, name);
  auto loc = local_buffer();
  auto loc_x = x.local_buffer();
  if (!loc->compatible(*loc_x))
    error(name + " incompatible local buffers");
  for (size_t i = 0; i < loc->size(); ++i)
    loc->at(i) *= loc_x->at(i);
}
void DistrArray::times(const DistrArray& x, const DistrArray& y) {
  auto name = std::string{"Array::times"};
  if (!compatible(x))
    error(name + " array x is incompatible");
  if (!compatible(y))
    error(name + " array y is incompatible");
  auto p = util::ScopeProfiler(m_prof, name);
  auto loc = local_buffer();
  auto loc_x = x.local_buffer();
  auto loc_y = x.local_buffer();
  if (!loc->compatible(*loc_x) || !loc->compatible(*loc_y))
    error(name + " incompatible local buffers");
  for (size_t i = 0; i < loc->size(); ++i)
    loc->at(i) = loc_x->at(i) * loc_y->at(i);
}

DistrArray::value_type DistrArray::dot(const DistrArray& x) const {
  auto name = std::string{"Array::dot"};
  if (!compatible(x))
    error(name + " array x is incompatible");
  if (empty() || x.empty())
    error(name + " calling dot on empty arrays");
  auto p = util::ScopeProfiler(m_prof, name);
  auto loc = local_buffer();
  auto loc_x = x.local_buffer();
  if (!loc->compatible(*loc_x))
    error(name + " incompatible local buffers");
  double a = std::inner_product(loc->begin(), loc->end(), loc_x->begin(), 0.);
  MPI_Allreduce(MPI_IN_PLACE, &a, 1, MPI_DOUBLE, MPI_SUM, m_communicator);
  return a;
}

void DistrArray::_divide(const DistrArray& x, const DistrArray& y, DistrArray::value_type shift, bool append,
                        bool negative) {
  auto name = std::string{"Array::divide"};
  if (!compatible(x))
    error(name + " array x is incompatible");
  if (empty() || x.empty())
    error(name + " calling dot on empty arrays");
  auto p = util::ScopeProfiler(m_prof, name);
  auto loc = local_buffer();
  auto loc_x = x.local_buffer();
  auto loc_y = y.local_buffer();
  if (!loc->compatible(*loc_x) || !loc->compatible(*loc_y))
    error(name + " incompatible local buffers");
  if (append) {
    if (negative)
      for (size_t i = 0; i < loc->size(); ++i)
        loc->at(i) -= loc_x->at(i) / (loc_y->at(i) + shift);
    else
      for (size_t i = 0; i < loc->size(); ++i)
        loc->at(i) += loc_x->at(i) / (loc_y->at(i) + shift);
  } else {
    if (negative)
      for (size_t i = 0; i < loc->size(); ++i)
        loc->at(i) = -loc_x->at(i) / (loc_y->at(i) + shift);
    else
      for (size_t i = 0; i < loc->size(); ++i)
        loc->at(i) = loc_x->at(i) / (loc_y->at(i) + shift);
  }
}

std::list<std::pair<size_t, DistrArray::value_type>> DistrArray::min_n(size_t n) const {
  return extrema<std::less<double>>(n);
}

std::list<std::pair<size_t, DistrArray::value_type>> DistrArray::max_n(size_t n) const {
  return extrema<std::greater<double>>(n);
}
std::list<std::pair<size_t, DistrArray::value_type>> DistrArray::min_abs_n(size_t n) const {
  return extrema<util::CompareAbs<double, std::less<>>>(n);
}
std::list<std::pair<size_t, DistrArray::value_type>> DistrArray::max_abs_n(size_t n) const {
  return extrema<util::CompareAbs<double, std::greater<>>>(n);
}
std::vector<size_t> DistrArray::min_loc_n(size_t n) const {
  auto min_list = min_abs_n(n);
  auto min_vec = std::vector<size_t>(n);
  std::transform(min_list.cbegin(), min_list.cend(), min_vec.begin(), [](const auto& p) { return p.first; });
  return min_vec;
}

template <class Compare>
std::list<std::pair<DistrArray::index_type, DistrArray::value_type>> DistrArray::extrema(int n) const {
  if (empty())
    return {};
  auto prof = util::ScopeProfiler(m_prof, "Array::extrema");
  auto buffer = local_buffer();
  auto length = buffer->size();
  auto nmin = length > n ? n : length;
  auto loc_extrema = std::list<std::pair<DistrArray::index_type, double>>();
  for (size_t i = 0; i < nmin; ++i)
    loc_extrema.emplace_back(buffer->lo + i, buffer->at(i));
  auto compare = Compare();
  auto compare_pair = [&compare](const auto& p1, const auto& p2) { return compare(p1.second, p2.second); };
  for (size_t i = nmin; i < length; ++i) {
    loc_extrema.emplace_back(buffer->lo + i, buffer->at(i));
    loc_extrema.sort(compare_pair);
    loc_extrema.pop_back();
  }
  auto indices_loc = std::vector<DistrArray::index_type>(n, size() + 1);
  auto indices_glob = std::vector<DistrArray::index_type>(n);
  auto values_loc = std::vector<double>(n);
  auto values_glob = std::vector<double>(n);
  size_t ind = 0;
  for (auto [i, v] : loc_extrema) {
    indices_loc[ind] = i;
    values_loc[ind] = v;
    ++ind;
  }
  MPI_Request requests[3];
  int comm_rank, comm_size;
  MPI_Comm_rank(m_communicator, &comm_rank);
  MPI_Comm_size(m_communicator, &comm_size);
  // root collects values, does the final sort and sends the result back
  if (comm_rank == 0) {
    auto ntot = n * comm_size;
    indices_loc.resize(ntot);
    values_loc.resize(ntot);
    auto ndummy = std::vector<int>(comm_size);
    auto d = int(n - nmin);
    MPI_Igather(&d, 1, MPI_INT, ndummy.data(), 1, MPI_INT, 0, m_communicator, &requests[0]);
    MPI_Igather(MPI_IN_PLACE, n, MPI_UNSIGNED_LONG, indices_loc.data(), n, MPI_UNSIGNED_LONG, 0, m_communicator,
                &requests[1]);
    MPI_Igather(MPI_IN_PLACE, n, MPI_DOUBLE, values_loc.data(), n, MPI_UNSIGNED_LONG, 0, m_communicator, &requests[2]);
    MPI_Waitall(3, requests, MPI_STATUSES_IGNORE);
    auto tot_dummy = std::accumulate(ndummy.cbegin(), ndummy.cend(), 0);
    if (tot_dummy != 0) {
      size_t shift = 0;
      for (size_t i = 0, ind = 0; i < comm_size; ++i) {
        for (size_t j = 0; j < n - ndummy[i]; ++j, ++ind) {
          indices_loc[ind] = indices_loc[ind + shift];
          values_loc[ind] = values_loc[ind + shift];
        }
        shift += ndummy[i];
      }
      indices_loc.resize(ntot - tot_dummy);
      values_loc.resize(ntot - tot_dummy);
    }
    std::vector<double> sort_permutation(indices_loc.size());
    std::iota(sort_permutation.begin(), sort_permutation.end(), 0);
    std::sort(
        sort_permutation.begin(), sort_permutation.end(),
        [&values_loc, &compare](const auto& i1, const auto& i2) { return compare(values_loc[i1], values_loc[i2]); });
    for (size_t i = 0; i < n; ++i) {
      auto j = sort_permutation[i];
      indices_glob[i] = indices_loc[j];
      values_glob[i] = values_loc[j];
    }
  } else {
    auto d = int(n - nmin);
    MPI_Igather(&d, 1, MPI_INT, nullptr, 1, MPI_INT, 0, m_communicator, &requests[0]);
    MPI_Igather(indices_loc.data(), n, MPI_UNSIGNED_LONG, nullptr, n, MPI_UNSIGNED_LONG, 0, m_communicator,
                &requests[1]);
    MPI_Igather(values_loc.data(), n, MPI_DOUBLE, nullptr, n, MPI_DOUBLE, 0, m_communicator, &requests[2]);
    MPI_Waitall(3, requests, MPI_STATUSES_IGNORE);
  }
  MPI_Ibcast(indices_glob.data(), n, MPI_UNSIGNED_LONG, 0, m_communicator, &requests[0]);
  MPI_Ibcast(values_glob.data(), n, MPI_DOUBLE, 0, m_communicator, &requests[1]);
  MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);
  auto map_extrema = std::list<std::pair<DistrArray::index_type, DistrArray::value_type>>();
  for (size_t i = 0; i < n; ++i)
    map_extrema.emplace_back(indices_glob[i], values_glob[i]);
  return map_extrema;
}
} // namespace molpro::gci::array