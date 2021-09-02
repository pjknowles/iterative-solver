#include "DistrArray.h"
#include "util/select.h"
#include "util/select_max_dot.h"
#include <algorithm>
#include <functional>
#include <iostream>
#include <numeric>
#include <molpro/Profiler.h>

namespace molpro::linalg::array {

DistrArray::DistrArray(size_t dimension, MPI_Comm commun) : m_dimension(dimension), m_communicator(commun) {}

void DistrArray::sync() const { MPI_Barrier(m_communicator); }

void DistrArray::error(const std::string& message) const {
  std::cerr << message << std::endl;
#ifdef HAVE_MPI_H
  MPI_Abort(m_communicator, 1);
#else
  throw std::runtime_error(message);
#endif
}

bool DistrArray::compatible(const DistrArray& other) const {
  bool result = (m_dimension == other.m_dimension);
#ifdef HAVE_MPI_H
  if (m_communicator == other.m_communicator)
    result &= true;
  else if (m_communicator == MPI_COMM_NULL || other.m_communicator == MPI_COMM_NULL)
    result &= false;
  else {
    int comp;
    MPI_Comm_compare(m_communicator, other.m_communicator, &comp);
    result &= (comp == MPI_IDENT || comp == MPI_CONGRUENT);
  }
#endif
  return result;
}

void DistrArray::zero() { fill(0); }

void DistrArray::fill(DistrArray::value_type val) {
  auto lb = local_buffer();
  for (auto& el : *lb)
    el = val;
}

void DistrArray::axpy(value_type a, const DistrArray& y) {
  auto prof = molpro::Profiler::single();
  prof->start("DistrArray::axpy");
  auto name = std::string{"Array::axpy"};
  if (!compatible(y))
    error(name + " incompatible arrays");
  if (a == 0)
    return;
  auto loc_x = local_buffer();
  auto loc_y = y.local_buffer();
  if (!loc_x->compatible(*loc_y))
    error(name + " incompatible local buffers");
  if (a == 1)
    for (size_t i = 0; i < loc_x->size(); ++i)
      (*loc_x)[i] += (*loc_y)[i];
  else if (a == -1)
    for (size_t i = 0; i < loc_x->size(); ++i)
      (*loc_x)[i] -= (*loc_y)[i];
  else
    for (size_t i = 0; i < loc_x->size(); ++i)
      (*loc_x)[i] += a * (*loc_y)[i];
  prof->stop();
}

void DistrArray::scal(DistrArray::value_type a) {
  auto x = local_buffer();
  for (auto& el : *x)
    el *= a;
}

void DistrArray::add(const DistrArray& y) { return axpy(1, y); }

void DistrArray::add(DistrArray::value_type a) {
  auto x = local_buffer();
  for (auto& el : *x)
    el += a;
}

void DistrArray::sub(const DistrArray& y) { return axpy(-1, y); }

void DistrArray::sub(DistrArray::value_type a) { return add(-a); }

void DistrArray::recip() {
  auto x = local_buffer();
  for (auto& el : *x)
    el = 1. / el;
}

void DistrArray::times(const DistrArray& y) {
  auto name = std::string{"Array::times"};
  if (!compatible(y))
    error(name + " incompatible arrays");
  auto loc_x = local_buffer();
  auto loc_y = y.local_buffer();
  if (!loc_x->compatible(*loc_y))
    error(name + " incompatible local buffers");
  for (size_t i = 0; i < loc_x->size(); ++i)
    (*loc_x)[i] *= (*loc_y)[i];
}

void DistrArray::times(const DistrArray& y, const DistrArray& z) {
  auto name = std::string{"Array::times"};
  if (!compatible(y))
    error(name + " array y is incompatible");
  if (!compatible(z))
    error(name + " array z is incompatible");
  auto loc_x = local_buffer();
  auto loc_y = y.local_buffer();
  auto loc_z = z.local_buffer();
  if (!loc_x->compatible(*loc_y) || !loc_x->compatible(*loc_z))
    error(name + " incompatible local buffers");
  for (size_t i = 0; i < loc_x->size(); ++i)
    (*loc_x)[i] = (*loc_y)[i] * (*loc_z)[i];
}

DistrArray::value_type DistrArray::dot(const DistrArray& y) const {
  auto prof = molpro::Profiler::single()->push("DistArray::dot()");
  auto name = std::string{"Array::dot"};
  if (!compatible(y))
    error(name + " array x is incompatible");
  auto loc_x = local_buffer();
  auto loc_y = y.local_buffer();
  if (!loc_x->compatible(*loc_y))
    error(name + " incompatible local buffers");
  auto a = std::inner_product(begin(*loc_x), end(*loc_x), begin(*loc_y), (value_type)0);
#ifdef HAVE_MPI_H
  MPI_Allreduce(MPI_IN_PLACE, &a, 1, MPI_DOUBLE, MPI_SUM, communicator());
#endif
  return a;
}

void DistrArray::_divide(const DistrArray& y, const DistrArray& z, DistrArray::value_type shift, bool append,
                         bool negative) {
  auto name = std::string{"Array::divide"};
  if (!compatible(y))
    error(name + " array y is incompatible");
  if (!compatible(z))
    error(name + " array z is incompatible");
  auto loc_x = local_buffer();
  auto loc_y = y.local_buffer();
  auto loc_z = z.local_buffer();
  if (!loc_x->compatible(*loc_y) || !loc_x->compatible(*loc_z))
    error(name + " incompatible local buffers");
  if (append) {
    if (negative)
      for (size_t i = 0; i < loc_x->size(); ++i)
        (*loc_x)[i] -= (*loc_y)[i] / ((*loc_z)[i] + shift);
    else
      for (size_t i = 0; i < loc_x->size(); ++i)
        (*loc_x)[i] += (*loc_y)[i] / ((*loc_z)[i] + shift);
  } else {
    if (negative)
      for (size_t i = 0; i < loc_x->size(); ++i)
        (*loc_x)[i] = -(*loc_y)[i] / ((*loc_z)[i] + shift);
    else
      for (size_t i = 0; i < loc_x->size(); ++i)
        (*loc_x)[i] = (*loc_y)[i] / ((*loc_z)[i] + shift);
  }
}

namespace util {
std::map<size_t, double> select_max_dot_broadcast(size_t n, std::map<size_t, double>& local_selection,
                                                  MPI_Comm communicator) {
#ifdef HAVE_MPI_H
  auto indices = std::vector<DistrArray::index_type>();
  auto values = std::vector<double>();
  indices.reserve(n);
  values.reserve(n);
  for (const auto& el : local_selection) {
    indices.push_back(el.first);
    values.push_back(el.second);
  }
  indices.resize(n);
  values.resize(n);
  int n_dummy = static_cast<int>(n) - local_selection.size();
  MPI_Request requests[3];
  int comm_rank, comm_size;
  MPI_Comm_rank(communicator, &comm_rank);
  MPI_Comm_size(communicator, &comm_size);
  if (comm_rank == 0) {
    auto n_tot = n * comm_size;
    auto n_dummy_elements = std::vector<int>(comm_size, 0);
    MPI_Igather(&n_dummy, 1, MPI_INT, &n_dummy_elements[0], 1, MPI_INT, 0, communicator, &requests[0]);
    indices.resize(n_tot);
    values.resize(n_tot);
    MPI_Igather(MPI_IN_PLACE, n, MPI_UNSIGNED_LONG, &indices[0], n, MPI_UNSIGNED_LONG, 0, communicator, &requests[1]);
    MPI_Igather(MPI_IN_PLACE, n, MPI_DOUBLE, &values[0], n, MPI_UNSIGNED_LONG, 0, communicator, &requests[2]);
    MPI_Waitall(3, requests, MPI_STATUSES_IGNORE);
    local_selection.clear();
    using pair_t = std::pair<double, size_t>;
    auto pq = std::priority_queue<pair_t, std::vector<pair_t>, std::greater<>>();
    for (size_t i = 0; i < n; ++i)
      pq.emplace(-std::numeric_limits<double>::max(), std::numeric_limits<DistrArray::index_type>::max());
    for (int i = 0, ii = 0; i < comm_size; ++i) {
      for (size_t j = 0; j < n - n_dummy_elements[i]; ++j, ++ii) {
        pq.emplace(values[ii], indices[ii]);
        pq.pop();
      }
      ii += n_dummy_elements[i];
    }
    indices.resize(n);
    values.resize(n);
    for (size_t i = 0; i < n; ++i) {
      values[i] = pq.top().first;
      indices[i] = pq.top().second;
      pq.pop();
    }
  } else {
    MPI_Igather(&n_dummy, 1, MPI_INT, nullptr, 1, MPI_INT, 0, communicator, &requests[0]);
    MPI_Igather(&indices[0], n, MPI_UNSIGNED_LONG, nullptr, n, MPI_UNSIGNED_LONG, 0, communicator, &requests[1]);
    MPI_Igather(&values[0], n, MPI_DOUBLE, nullptr, n, MPI_DOUBLE, 0, communicator, &requests[2]);
    MPI_Waitall(3, requests, MPI_STATUSES_IGNORE);
  }
  MPI_Ibcast(&indices[0], n, MPI_UNSIGNED_LONG, 0, communicator, &requests[0]);
  MPI_Ibcast(&values[0], n, MPI_DOUBLE, 0, communicator, &requests[1]);
  MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);
  local_selection.clear();
  for (size_t i = 0; i < n; ++i)
    local_selection.emplace(indices[i], values[i]);
#endif
  return local_selection;
}
} // namespace util

std::map<size_t, DistrArray::value_type> DistrArray::select_max_dot(size_t n, const DistrArray& y) const {
  if (!compatible(y))
    error("DistrArray::select_max_dot: incompatible arrays");
  if (n > size() || n > y.size())
    error("DistrArray::select_max_dot: n is too large");
  auto xbuf = local_buffer();
  auto ybuf = y.local_buffer();
  auto local_selection =
      util::select_max_dot<LocalBuffer, LocalBuffer, value_type, value_type>(std::min(n, xbuf->size()), *xbuf, *ybuf);
  auto shifted_local_selection = decltype(local_selection)();
  for (auto& el : local_selection)
    shifted_local_selection.emplace(xbuf->start() + el.first, el.second);
  return util::select_max_dot_broadcast(n, shifted_local_selection, communicator());
}

std::map<size_t, DistrArray::value_type> DistrArray::select_max_dot(size_t n, const DistrArray::SparseArray& y) const {
  auto name = std::string("DistrArray::select_max_dot:");
  if (size() < y.rbegin()->first + 1)
    error(name + " sparse array x is too large");
  if (n > size() || n > y.size())
    error(" n is too large");
  auto xbuf = local_buffer();
  auto local_selection = util::select_max_dot_iter_sparse<LocalBuffer, SparseArray, value_type, value_type>(
      std::min(n, xbuf->size()), *xbuf, y);
  auto shifted_local_selection = decltype(local_selection)();
  for (auto& el : local_selection)
    shifted_local_selection.emplace(xbuf->start() + el.first, el.second);
  return util::select_max_dot_broadcast(n, shifted_local_selection, communicator());
}

std::map<size_t, DistrArray::value_type> DistrArray::select(size_t n, bool max, bool ignore_sign) const {
  if (n > size())
    error("DistrArray::select: n is too large");
  auto xbuf = local_buffer();
  auto local_selection = util::select<LocalBuffer, value_type>(std::min(n, xbuf->size()), *xbuf, max, ignore_sign);
  auto shifted_local_selection = decltype(local_selection)();
  for (const auto& el : local_selection)
    shifted_local_selection.emplace(xbuf->start() + el.first, max ? el.second : -el.second);
  std::map<size_t, double> result = util::select_max_dot_broadcast(n, shifted_local_selection, communicator());
  if (not max)
    for (auto& el : result)
      el.second = -el.second;
  return result;
}

namespace util {
template <class Compare>
std::list<std::pair<DistrArray::index_type, DistrArray::value_type>> extrema(const DistrArray& x, int n) {
  if (x.size() == 0)
    return {};
  auto buffer = x.local_buffer();
  auto length = buffer->size();
  auto nmin = length > size_t(n) ? size_t(n) : length;
  auto loc_extrema = std::list<std::pair<DistrArray::index_type, double>>();
  for (size_t i = 0; i < nmin; ++i)
    loc_extrema.emplace_back(buffer->start() + i, (*buffer)[i]);
  auto compare = Compare();
  auto compare_pair = [&compare](const auto& p1, const auto& p2) { return compare(p1.second, p2.second); };
  for (size_t i = nmin; i < length; ++i) {
    loc_extrema.emplace_back(buffer->start() + i, (*buffer)[i]);
    loc_extrema.sort(compare_pair);
    loc_extrema.pop_back();
  }
  auto indices_loc = std::vector<DistrArray::index_type>(n, x.size() + 1);
  auto indices_glob = std::vector<DistrArray::index_type>(n);
  auto values_loc = std::vector<double>(n);
  auto values_glob = std::vector<double>(n);
  size_t ind = 0;
  for (auto& it : loc_extrema) {
    indices_loc[ind] = it.first;
    values_loc[ind] = it.second;
    ++ind;
  }
#ifdef HAVE_MPI_H
  MPI_Request requests[3];
  int comm_rank, comm_size;
  MPI_Comm_rank(x.communicator(), &comm_rank);
  MPI_Comm_size(x.communicator(), &comm_size);
  // root collects values, does the final sort and sends the result back
  if (comm_rank == 0) {
    auto ntot = n * comm_size;
    indices_loc.resize(ntot);
    values_loc.resize(ntot);
    auto ndummy = std::vector<int>(comm_size);
    auto d = int(n - nmin);
    MPI_Igather(&d, 1, MPI_INT, ndummy.data(), 1, MPI_INT, 0, x.communicator(), &requests[0]);
    MPI_Igather(MPI_IN_PLACE, n, MPI_UNSIGNED_LONG, indices_loc.data(), n, MPI_UNSIGNED_LONG, 0, x.communicator(),
                &requests[1]);
    MPI_Igather(MPI_IN_PLACE, n, MPI_DOUBLE, values_loc.data(), n, MPI_UNSIGNED_LONG, 0, x.communicator(),
                &requests[2]);
    MPI_Waitall(3, requests, MPI_STATUSES_IGNORE);
    auto tot_dummy = std::accumulate(begin(ndummy), end(ndummy), 0);
    if (tot_dummy != 0) {
      size_t shift = 0;
      for (int i = 0, ind = 0; i < comm_size; ++i) {
        for (int j = 0; j < n - ndummy[i]; ++j, ++ind) {
          indices_loc[ind] = indices_loc[ind + shift];
          values_loc[ind] = values_loc[ind + shift];
        }
        shift += ndummy[i];
      }
      indices_loc.resize(ntot - tot_dummy);
      values_loc.resize(ntot - tot_dummy);
    }
#endif
    std::vector<unsigned int> sort_permutation(indices_loc.size());
    std::iota(begin(sort_permutation), end(sort_permutation), (unsigned int)0);
    std::sort(begin(sort_permutation), end(sort_permutation), [&values_loc, &compare](const auto& i1, const auto& i2) {
      return compare(values_loc[i1], values_loc[i2]);
    });
    for (int i = 0; i < n; ++i) {
      auto j = sort_permutation[i];
      indices_glob[i] = indices_loc[j];
      values_glob[i] = values_loc[j];
    }
#ifdef HAVE_MPI_H
  } else {
    auto d = int(n - nmin);
    MPI_Igather(&d, 1, MPI_INT, nullptr, 1, MPI_INT, 0, x.communicator(), &requests[0]);
    MPI_Igather(indices_loc.data(), n, MPI_UNSIGNED_LONG, nullptr, n, MPI_UNSIGNED_LONG, 0, x.communicator(),
                &requests[1]);
    MPI_Igather(values_loc.data(), n, MPI_DOUBLE, nullptr, n, MPI_DOUBLE, 0, x.communicator(), &requests[2]);
    MPI_Waitall(3, requests, MPI_STATUSES_IGNORE);
  }
  MPI_Ibcast(indices_glob.data(), n, MPI_UNSIGNED_LONG, 0, x.communicator(), &requests[0]);
  MPI_Ibcast(values_glob.data(), n, MPI_DOUBLE, 0, x.communicator(), &requests[1]);
  MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);
#endif
  auto map_extrema = std::list<std::pair<DistrArray::index_type, double>>();
  for (int i = 0; i < n; ++i)
    map_extrema.emplace_back(indices_glob[i], values_glob[i]);
  return map_extrema;
}
} // namespace util

std::list<std::pair<DistrArray::index_type, DistrArray::value_type>> DistrArray::min_n(int n) const {
  return util::extrema<std::less<DistrArray::value_type>>(*this, n);
}

std::list<std::pair<DistrArray::index_type, DistrArray::value_type>> DistrArray::max_n(int n) const {
  return util::extrema<std::greater<DistrArray::value_type>>(*this, n);
}

std::list<std::pair<DistrArray::index_type, DistrArray::value_type>> DistrArray::min_abs_n(int n) const {
  return util::extrema<util::CompareAbs<DistrArray::value_type, std::less<>>>(*this, n);
}

std::list<std::pair<DistrArray::index_type, DistrArray::value_type>> DistrArray::max_abs_n(int n) const {
  return util::extrema<util::CompareAbs<DistrArray::value_type, std::greater<>>>(*this, n);
}

std::vector<DistrArray::index_type> DistrArray::min_loc_n(int n) const {
  auto min_list = min_abs_n(n);
  auto min_vec = std::vector<index_type>(n);
  std::transform(begin(min_list), end(min_list), begin(min_vec), [](const auto& p) { return p.first; });
  return min_vec;
}

void DistrArray::copy(const DistrArray& y) {
  auto name = std::string{"Array::copy"};
  if (!compatible(y))
    error(name + " incompatible arrays");
  auto loc_x = local_buffer();
  auto loc_y = y.local_buffer();
  if (!loc_x->compatible(*loc_y))
    error(name + " incompatible local buffers");
  for (size_t i = 0; i < loc_x->size(); ++i)
    (*loc_x)[i] = (*loc_y)[i];
}

void DistrArray::copy_patch(const DistrArray& y, DistrArray::index_type start, DistrArray::index_type end) {
  auto name = std::string{"Array::copy_patch"};
  if (!compatible(y))
    error(name + " incompatible arrays");
  auto loc_x = local_buffer();
  auto loc_y = y.local_buffer();
  if (!loc_x->compatible(*loc_y))
    error(name + " incompatible local buffers");
  if (start > end)
    return;
  auto s = start <= loc_x->start() ? 0 : start - loc_x->start();
  auto e = end - start + 1 >= loc_x->size() ? loc_x->size() : end - start + 1;
  for (auto i = s; i < e; ++i)
    (*loc_x)[i] = (*loc_y)[i];
}

DistrArray::value_type DistrArray::dot(const SparseArray& y) const {
  auto name = std::string{"Array::dot SparseArray "};
  if (size() < y.rbegin()->first + 1)
    error(name + " sparse array x is incompatible");
  auto loc_x = local_buffer();
  double res = 0;
  if (loc_x->size() > 0) {
    index_type i;
    value_type v;
    for (auto it = y.lower_bound(loc_x->start()); it != y.upper_bound(loc_x->start() + loc_x->size() - 1); ++it) {
      std::tie(i, v) = *it;
      res += (*loc_x)[i - loc_x->start()] * v;
    }
  }
#ifdef HAVE_MPI_H
  MPI_Allreduce(MPI_IN_PLACE, &res, 1, MPI_DOUBLE, MPI_SUM, communicator());
#endif
  return res;
}

void DistrArray::axpy(value_type a, const SparseArray& y) {
  auto name = std::string{"Array::axpy SparseArray"};
  if (a == 0)
    return;
  if (size() < y.rbegin()->first + 1)
    error(name + " sparse array x is incompatible");
  auto loc_x = local_buffer();
  if (loc_x->size() > 0) {
    index_type i;
    value_type v;
    if (a == 1)
      for (auto it = y.lower_bound(loc_x->start()); it != y.upper_bound(loc_x->start() + loc_x->size() - 1); ++it) {
        std::tie(i, v) = *it;
        (*loc_x)[i - loc_x->start()] += v;
      }
    else if (a == -1)
      for (auto it = y.lower_bound(loc_x->start()); it != y.upper_bound(loc_x->start() + loc_x->size() - 1); ++it) {
        std::tie(i, v) = *it;
        (*loc_x)[i - loc_x->start()] -= v;
      }
    else
      for (auto it = y.lower_bound(loc_x->start()); it != y.upper_bound(loc_x->start() + loc_x->size() - 1); ++it) {
        std::tie(i, v) = *it;
        (*loc_x)[i - loc_x->start()] += a * v;
      }
  }
}

} // namespace molpro::linalg::array