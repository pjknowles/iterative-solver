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

bool DistrArray::compatible(const DistrArray& other) const {
  int comp;
  MPI_Comm_compare(m_communicator, other.m_communicator, &comp);
  return (m_dimension == other.m_dimension) && (MPI_IDENT || MPI_CONGRUENT);
}

bool DistrArray::empty() const { return true; }

void DistrArray::zero() { fill(*this, 0); }

void fill(DistrArray& x, DistrArray::value_type val) {
  auto p = util::ScopeProfiler(x.m_prof, "Array::fill");
  for (auto& el : *x.local_buffer())
    el = val;
}

void axpy(DistrArray& x, DistrArray::value_type a, const DistrArray& y) {
  auto name = std::string{"Array::axpy"};
  if (!x.compatible(y))
    x.error(name + " incompatible arrays");
  if (x.empty() || y.empty())
    x.error(name + " cannot use empty arrays");
  if (a == 0)
    return;
  auto p = util::ScopeProfiler(x.m_prof, name);
  auto loc_x = x.local_buffer();
  auto loc_y = y.local_buffer();
  if (!loc_x->compatible(*loc_y))
    x.error(name + " incompatible local buffers");
  if (a == 1)
    for (size_t i = 0; i < loc_x->size(); ++i)
      loc_x->at(i) += loc_y->at(i);
  else if (a == -1)
    for (size_t i = 0; i < loc_x->size(); ++i)
      loc_x->at(i) -= loc_y->at(i);
  else
    for (size_t i = 0; i < loc_x->size(); ++i)
      loc_x->at(i) += a * loc_y->at(i);
}

void scal(DistrArray& x, DistrArray::value_type a) {
  auto p = util::ScopeProfiler(x.m_prof, "Array::scal");
  for (auto& el : *x.local_buffer())
    el *= a;
}

void add(DistrArray& x, const DistrArray& y) { axpy(x, 1, y); }

void add(DistrArray& x, DistrArray::value_type a);

void add(DistrArray& x, DistrArray::value_type a) {
  auto p = util::ScopeProfiler(x.m_prof, "Array::add");
  for (auto& el : *x.local_buffer())
    el += a;
}

void sub(DistrArray& x, const DistrArray& y) { axpy(x, -1, y); }

void sub(DistrArray& x, DistrArray::value_type a) { add(x, -a); }

void recip(DistrArray& x) {
  auto p = util::ScopeProfiler(x.m_prof, "Array::recip");
  for (auto& el : *x.local_buffer())
    el = 1. / el;
}

void times(DistrArray& x, const DistrArray& y) {
  auto name = std::string{"Array::times"};
  if (!x.compatible(y))
    x.error(name + " incompatible arrays");
  if (x.empty() || y.empty())
    x.error(name + " cannot use empty arrays");
  auto p = util::ScopeProfiler(x.m_prof, name);
  auto loc_x = x.local_buffer();
  auto loc_y = y.local_buffer();
  if (!loc_x->compatible(*loc_y))
    x.error(name + " incompatible local buffers");
  for (size_t i = 0; i < loc_x->size(); ++i)
    loc_x->at(i) *= loc_y->at(i);
}

void times(DistrArray& x, const DistrArray& y, const DistrArray& z) {
  auto name = std::string{"Array::times"};
  if (!x.compatible(y))
    x.error(name + " array y is incompatible");
  if (!x.compatible(z))
    x.error(name + " array z is incompatible");
  if (x.empty() || y.empty() || z.empty())
    x.error(name + " cannot use empty arrays");
  auto p = util::ScopeProfiler(x.m_prof, name);
  auto loc_x = x.local_buffer();
  auto loc_y = y.local_buffer();
  auto loc_z = z.local_buffer();
  if (!loc_x->compatible(*loc_y) || !loc_x->compatible(*loc_z))
    x.error(name + " incompatible local buffers");
  for (size_t i = 0; i < loc_x->size(); ++i)
    loc_x->at(i) = loc_y->at(i) * loc_z->at(i);
}

DistrArray::value_type dot(const DistrArray& x, const DistrArray& y) {
  auto name = std::string{"Array::dot"};
  if (!x.compatible(y))
    x.error(name + " array x is incompatible");
  if (x.empty() || y.empty())
    x.error(name + " calling dot on empty arrays");
  auto p = util::ScopeProfiler(x.m_prof, name);
  auto loc_x = x.local_buffer();
  auto loc_y = y.local_buffer();
  if (!loc_x->compatible(*loc_y))
    x.error(name + " incompatible local buffers");
  DistrArray::value_type a = std::inner_product(loc_x->begin(), loc_x->end(), loc_y->begin(), 0.);
  MPI_Allreduce(MPI_IN_PLACE, &a, 1, MPI_DOUBLE, MPI_SUM, x.communicator());
  return a;
}

void divide(DistrArray& x, const DistrArray& y, const DistrArray& z, DistrArray::value_type shift, bool append,
            bool negative) {
  auto name = std::string{"Array::divide"};
  if (!x.compatible(y))
    x.error(name + " array y is incompatible");
  if (!x.compatible(z))
    x.error(name + " array z is incompatible");
  if (x.empty() || y.empty() || z.empty())
    x.error(name + " calling divide with an empty array");
  auto p = util::ScopeProfiler(x.m_prof, name);
  auto loc_x = x.local_buffer();
  auto loc_y = y.local_buffer();
  auto loc_z = z.local_buffer();
  if (!loc_x->compatible(*loc_y) || !loc_x->compatible(*loc_z))
    x.error(name + " incompatible local buffers");
  if (append) {
    if (negative)
      for (size_t i = 0; i < loc_x->size(); ++i)
        loc_x->at(i) -= loc_y->at(i) / (loc_z->at(i) + shift);
    else
      for (size_t i = 0; i < loc_x->size(); ++i)
        loc_x->at(i) += loc_y->at(i) / (loc_z->at(i) + shift);
  } else {
    if (negative)
      for (size_t i = 0; i < loc_x->size(); ++i)
        loc_x->at(i) = -loc_y->at(i) / (loc_z->at(i) + shift);
    else
      for (size_t i = 0; i < loc_x->size(); ++i)
        loc_x->at(i) = loc_y->at(i) / (loc_z->at(i) + shift);
  }
}

namespace util {
template <class Compare> std::list<std::pair<unsigned long, double>> extrema(const DistrArray& x, int n) {
  if (x.empty())
    return {};
  auto prof = util::ScopeProfiler(x.m_prof, "Array::extrema");
  auto buffer = x.local_buffer();
  auto length = buffer->size();
  auto nmin = length > n ? n : length;
  auto loc_extrema = std::list<std::pair<unsigned long, double>>();
  for (size_t i = 0; i < nmin; ++i)
    loc_extrema.emplace_back(buffer->lo + i, buffer->at(i));
  auto compare = Compare();
  auto compare_pair = [&compare](const auto& p1, const auto& p2) { return compare(p1.second, p2.second); };
  for (size_t i = nmin; i < length; ++i) {
    loc_extrema.emplace_back(buffer->lo + i, buffer->at(i));
    loc_extrema.sort(compare_pair);
    loc_extrema.pop_back();
  }
  auto indices_loc = std::vector<unsigned long>(n, x.size() + 1);
  auto indices_glob = std::vector<unsigned long>(n);
  auto values_loc = std::vector<unsigned long>(n);
  auto values_glob = std::vector<unsigned long>(n);
  size_t ind = 0;
  for (auto [i, v] : loc_extrema) {
    indices_loc[ind] = i;
    values_loc[ind] = v;
    ++ind;
  }
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
    MPI_Igather(&d, 1, MPI_INT, nullptr, 1, MPI_INT, 0, x.communicator(), &requests[0]);
    MPI_Igather(indices_loc.data(), n, MPI_UNSIGNED_LONG, nullptr, n, MPI_UNSIGNED_LONG, 0, x.communicator(),
                &requests[1]);
    MPI_Igather(values_loc.data(), n, MPI_DOUBLE, nullptr, n, MPI_DOUBLE, 0, x.communicator(), &requests[2]);
    MPI_Waitall(3, requests, MPI_STATUSES_IGNORE);
  }
  MPI_Ibcast(indices_glob.data(), n, MPI_UNSIGNED_LONG, 0, x.communicator(), &requests[0]);
  MPI_Ibcast(values_glob.data(), n, MPI_DOUBLE, 0, x.communicator(), &requests[1]);
  MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);
  auto map_extrema = std::list<std::pair<unsigned long, double>>();
  for (size_t i = 0; i < n; ++i)
    map_extrema.emplace_back(indices_glob[i], values_glob[i]);
  return map_extrema;
}
} // namespace util

std::list<std::pair<DistrArray::index_type, DistrArray::value_type>> min_n(const DistrArray& x, int n) {
  if (x.empty())
    return {};
  return util::extrema<std::less<DistrArray::value_type>>(x, n);
}

std::list<std::pair<DistrArray::index_type, DistrArray::value_type>> max_n(const DistrArray& x, int n) {
  if (x.empty())
    return {};
  return util::extrema<std::greater<DistrArray::value_type>>(x, n);
}

std::list<std::pair<DistrArray::index_type, DistrArray::value_type>> min_abs_n(const DistrArray& x, int n) {
  if (x.empty())
    return {};
  return util::extrema<util::CompareAbs<DistrArray::value_type, std::less<>>>(x, n);
}

std::list<std::pair<DistrArray::index_type, DistrArray::value_type>> max_abs_n(const DistrArray& x, int n) {
  if (x.empty())
    return {};
  return util::extrema<util::CompareAbs<DistrArray::value_type, std::greater<>>>(x, n);
}

std::vector<DistrArray::index_type> min_loc_n(const DistrArray& x, int n) {
  if (x.empty())
    return {};
  auto min_list = min_abs_n(x, n);
  auto min_vec = std::vector<DistrArray::index_type>(n);
  std::transform(min_list.cbegin(), min_list.cend(), min_vec.begin(), [](const auto& p) { return p.first; });
  return min_vec;
}

void copy(DistrArray& x, const DistrArray& y) {
  auto name = std::string{"Array::copy"};
  if (!x.compatible(y))
    x.error(name + " incompatible arrays");
  if (x.empty() != y.empty())
    x.error(name + " one of the arrays is empty");
  auto p = util::ScopeProfiler(x.m_prof, name);
  auto loc_x = x.local_buffer();
  auto loc_y = y.local_buffer();
  if (!loc_x->compatible(*loc_y))
    x.error(name + " incompatible local buffers");
  for (size_t i = 0; i < loc_x->size(); ++i)
    loc_x->at(i) = loc_y->at(i);
}
void copy_patch(DistrArray& x, const DistrArray& y, DistrArray::index_type start, DistrArray::index_type end) {
  auto name = std::string{"Array::copy_patch"};
  if (!x.compatible(y))
    x.error(name + " incompatible arrays");
  if (x.empty() != y.empty())
    x.error(name + " one of the arrays is empty");
  auto p = util::ScopeProfiler(x.m_prof, name);
  auto loc_x = x.local_buffer();
  auto loc_y = y.local_buffer();
  if (!loc_x->compatible(*loc_y))
    x.error(name + " incompatible local buffers");
  if (start > end)
    return;
  auto s = start <= loc_x->lo ? 0 : start - loc_x->lo;
  auto e = end - start + 1 >= loc_x->size() ? loc_x->size() : end - start + 1;
  for (auto i = s; i < e; ++i)
    loc_x->at(i) = loc_y->at(i);
}
} // namespace molpro::gci::array