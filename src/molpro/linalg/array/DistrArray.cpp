#include <molpro/linalg/scalar_traits.h>
#include "DistrArray.h"
#include "util/gather_all.h"
#include "util/select.h"
#include "util/select_max_dot.h"
#include <algorithm>
#include <functional>
#include <iostream>
#include <numeric>
#include <molpro/Profiler.h>

namespace molpro::linalg::array {

namespace {
using value_type = DistrArray::value_type;
using index_type = DistrArray::index_type;

//! One operand of an elementwise operation
struct Operand {
  //! A distributed array, of which the local section is operated on
  Operand(const DistrArray& array, bool read, bool write) : array(&array), read(read), write(write) {}
  //! An array held in full in contiguous memory, as many elements as the distributed operands, and addressed in place
  Operand(const value_type* memory, bool read, bool write)
      : memory(const_cast<value_type*>(memory)), read(read), write(write) {}
  const DistrArray* array = nullptr;
  value_type* memory = nullptr;
  bool read;  //!< whether the operation reads the operand
  bool write; //!< whether the operation writes the operand
};

//! The section [lo, hi) of the array held by this process
std::pair<index_type, index_type> local_range(const DistrArray& a) {
  int rank = 0;
  if (a.communicator() == molpro::mpi::comm_global())
    rank = molpro::mpi::rank_global();
#ifdef HAVE_MPI_H
  else
    MPI_Comm_rank(a.communicator(), &rank);
#endif
  return a.distribution().range(rank);
}

/*!
 * @brief Applies f(start, n, p) to the local sections of the operands, where p[k] addresses elements [start, start+n)
 * of operand k.
 *
 * The local section is that of the distributed operands, of which there must be at least one. Operands held in
 * memory, i.e. distributed arrays with a local buffer and contiguous arrays, are addressed in place, and if none of the
 * operands is held on disk, f is called once for the whole local section. Operands on disk, which have no local buffer,
 * are paged through buffers of disk_page_size() elements: each page is read before f is called if the operation reads
 * the operand, and written back afterwards if it writes it.
 */
void for_each_local_chunk(std::initializer_list<Operand> operands,
                          const std::function<void(index_type start, size_t n, value_type* const* p)>& f) {
  const std::vector<Operand> ops(operands);
  const auto& distributed = *std::find_if(ops.begin(), ops.end(), [](const auto& op) { return op.array; })->array;
  const auto [lo, hi] = local_range(distributed);
  // the local buffers of the operands in memory (DistrArray::LocalBuffer itself is not accessible here)
  std::vector<decltype(std::declval<DistrArray&>().local_buffer())> buffers;
  std::vector<decltype(std::declval<const DistrArray&>().local_buffer())> const_buffers;
  std::vector<value_type*> base(ops.size(), nullptr);
  std::vector<bool> on_disk(ops.size(), false);
  size_t page = 0;
  for (size_t k = 0; k < ops.size(); ++k) {
    if (not ops[k].array) {
      base[k] = ops[k].memory + lo; // an empty contiguous array may have no memory, but then lo == 0
      continue;
    }
    const auto& a = *ops[k].array;
    if (local_range(a) != std::make_pair(lo, hi))
      distributed.error("DistrArray: operands are distributed differently");
    if (a.disk_page_size() > 0) {
      on_disk[k] = true;
      page = page == 0 ? a.disk_page_size() : std::min(page, a.disk_page_size());
    } else if (ops[k].write) {
      buffers.push_back(const_cast<DistrArray&>(a).local_buffer());
      base[k] = buffers.back()->data();
    } else {
      const_buffers.push_back(a.local_buffer());
      base[k] = const_cast<value_type*>(const_buffers.back()->data());
    }
  }
  if (page == 0) {
    f(lo, hi - lo, base.data());
    return;
  }
  std::vector<std::vector<value_type>> pages(ops.size());
  std::vector<value_type*> p(ops.size());
  for (index_type start = lo; start < hi; start += page) {
    const size_t n = std::min<size_t>(page, hi - start);
    for (size_t k = 0; k < ops.size(); ++k) {
      if (on_disk[k]) {
        pages[k].resize(n);
        if (ops[k].read)
          ops[k].array->get(start, start + n, pages[k].data());
        p[k] = pages[k].data();
      } else
        p[k] = base[k] + (start - lo);
    }
    f(start, n, p.data());
    for (size_t k = 0; k < ops.size(); ++k)
      if (on_disk[k] and ops[k].write)
        const_cast<DistrArray&>(*ops[k].array).put(start, start + n, pages[k].data());
  }
}

//! Keeps the n entries of selection with the largest values
void keep_largest(std::map<size_t, value_type>& selection, size_t n) {
  if (selection.size() <= n)
    return;
  auto entries = std::vector<std::pair<size_t, value_type>>(selection.begin(), selection.end());
  std::nth_element(entries.begin(), entries.begin() + n, entries.end(),
                   [](const auto& a, const auto& b) { return a.second > b.second; });
  selection = std::map<size_t, value_type>(entries.begin(), entries.begin() + n);
}

//! The hermitian inner product <x|y> of two operands, reduced over the communicator
value_type dot(const Operand& x, const Operand& y, MPI_Comm communicator) {
  value_type a = 0;
  for_each_local_chunk({x, y}, [&a](index_type, size_t n, value_type* const* p) {
    a = std::inner_product(p[0], p[0] + n, p[1], a, std::plus<value_type>{}, [](const auto& elx, const auto& ely) {
      return molpro::linalg::conjugate(value_type(elx)) * value_type(ely);
    });
  });
#ifdef HAVE_MPI_H
  MPI_Allreduce(MPI_IN_PLACE, &a, 1, MPI_DOUBLE, MPI_SUM, communicator);
#endif
  return a;
}

//! x[:] += a * y[:] over the local section; x is the first operand
void axpy(value_type a, const Operand& x, const Operand& y) {
  if (a == 0)
    return;
  for_each_local_chunk({x, y}, [a](index_type, size_t n, value_type* const* p) {
    auto xx = p[0];
    const auto yy = p[1];
    if (a == 1)
      for (size_t i = 0; i < n; ++i)
        xx[i] += yy[i];
    else if (a == -1)
      for (size_t i = 0; i < n; ++i)
        xx[i] -= yy[i];
    else
      for (size_t i = 0; i < n; ++i)
        xx[i] += a * yy[i];
  });
}

//! x[:] = y[:] over the local section; x is the first operand
void copy(const Operand& x, const Operand& y) {
  for_each_local_chunk({x, y}, [](index_type, size_t n, value_type* const* p) { std::copy_n(p[1], n, p[0]); });
}

//! The n indices with the largest |x[i]*y[i]|, selected over all processes
std::map<size_t, value_type> select_max_dot(size_t n, const Operand& x, const Operand& y, MPI_Comm communicator) {
  auto shifted_local_selection = std::map<size_t, value_type>();
  for_each_local_chunk({x, y}, [n, &shifted_local_selection](index_type start, size_t m, value_type* const* p) {
    const auto xs = Span<value_type>(p[0], m);
    const auto ys = Span<value_type>(p[1], m);
    for (const auto& el :
         util::select_max_dot<Span<value_type>, Span<value_type>, value_type, value_type>(std::min(n, m), xs, ys))
      shifted_local_selection.emplace(start + el.first, el.second);
    keep_largest(shifted_local_selection, n);
  });
  return util::select_max_dot_broadcast(n, shifted_local_selection, communicator);
}
} // namespace

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

DistrArray::value_type DistrArray::operator[](index_type ind) const {
  return at(ind);
}

util::ValueProxy<DistrArray> DistrArray::operator[](index_type ind) {
  return {*this, std::move(ind)};
}

void DistrArray::zero() { fill(0); }

void DistrArray::fill(DistrArray::value_type val) {
  for_each_local_chunk({{*this, false, true}}, [val](index_type, size_t n, value_type* const* p) {
    std::fill_n(p[0], n, val);
  });
}

void DistrArray::axpy(value_type a, const DistrArray& y) {
  auto prof = molpro::Profiler::single();
  prof->start("DistrArray::axpy");
  auto name = std::string{"Array::axpy"};
  if (!compatible(y))
    error(name + " incompatible arrays");
  array::axpy(a, {*this, true, true}, {y, true, false});
  prof->stop();
}

void DistrArray::scal(DistrArray::value_type a) {
  for_each_local_chunk({{*this, true, true}}, [a](index_type, size_t n, value_type* const* p) {
    for (size_t i = 0; i < n; ++i)
      p[0][i] *= a;
  });
}

void DistrArray::add(const DistrArray& y) { return axpy(1, y); }

void DistrArray::add(DistrArray::value_type a) {
  for_each_local_chunk({{*this, true, true}}, [a](index_type, size_t n, value_type* const* p) {
    for (size_t i = 0; i < n; ++i)
      p[0][i] += a;
  });
}

void DistrArray::sub(const DistrArray& y) { return axpy(-1, y); }

void DistrArray::sub(DistrArray::value_type a) { return add(-a); }

void DistrArray::recip() {
  for_each_local_chunk({{*this, true, true}}, [](index_type, size_t n, value_type* const* p) {
    for (size_t i = 0; i < n; ++i)
      p[0][i] = 1. / p[0][i];
  });
}

void DistrArray::times(const DistrArray& y) {
  auto name = std::string{"Array::times"};
  if (!compatible(y))
    error(name + " incompatible arrays");
  for_each_local_chunk({{*this, true, true}, {y, true, false}}, [](index_type, size_t n, value_type* const* p) {
    for (size_t i = 0; i < n; ++i)
      p[0][i] *= p[1][i];
  });
}

void DistrArray::times(const DistrArray& y, const DistrArray& z) {
  auto name = std::string{"Array::times"};
  if (!compatible(y))
    error(name + " array y is incompatible");
  if (!compatible(z))
    error(name + " array z is incompatible");
  for_each_local_chunk({{*this, false, true}, {y, true, false}, {z, true, false}},
                       [](index_type, size_t n, value_type* const* p) {
                         for (size_t i = 0; i < n; ++i)
                           p[0][i] = p[1][i] * p[2][i];
                       });
}

DistrArray::value_type DistrArray::dot(const DistrArray& y) const {
  auto prof = molpro::Profiler::single()->push("DistrArray::dot()");
  auto name = std::string{"Array::dot"};
  if (!compatible(y))
    error(name + " array x is incompatible");
  return array::dot({*this, true, false}, {y, true, false}, communicator());
}

void DistrArray::_divide(const DistrArray& y, const DistrArray& z, DistrArray::value_type shift, bool append,
                         bool negative) {
  auto name = std::string{"Array::divide"};
  if (!compatible(y))
    error(name + " array y is incompatible");
  if (!compatible(z))
    error(name + " array z is incompatible");
  for_each_local_chunk({{*this, append, true}, {y, true, false}, {z, true, false}},
                       [shift, append, negative](index_type, size_t n, value_type* const* p) {
                         auto x = p[0];
                         const auto yy = p[1];
                         const auto zz = p[2];
                         if (append) {
                           if (negative)
                             for (size_t i = 0; i < n; ++i)
                               x[i] -= yy[i] / (zz[i] + shift);
                           else
                             for (size_t i = 0; i < n; ++i)
                               x[i] += yy[i] / (zz[i] + shift);
                         } else {
                           if (negative)
                             for (size_t i = 0; i < n; ++i)
                               x[i] = -yy[i] / (zz[i] + shift);
                           else
                             for (size_t i = 0; i < n; ++i)
                               x[i] = yy[i] / (zz[i] + shift);
                         }
                       });
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
    MPI_Igather(MPI_IN_PLACE, n, mpi_index_type(), &indices[0], n, mpi_index_type(), 0, communicator, &requests[1]);
    MPI_Igather(MPI_IN_PLACE, n, MPI_DOUBLE, &values[0], n, MPI_DOUBLE, 0, communicator, &requests[2]);
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
    MPI_Igather(&indices[0], n, mpi_index_type(), nullptr, n, mpi_index_type(), 0, communicator, &requests[1]);
    MPI_Igather(&values[0], n, MPI_DOUBLE, nullptr, n, MPI_DOUBLE, 0, communicator, &requests[2]);
    MPI_Waitall(3, requests, MPI_STATUSES_IGNORE);
  }
  MPI_Ibcast(&indices[0], n, mpi_index_type(), 0, communicator, &requests[0]);
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
  return array::select_max_dot(n, {*this, true, false}, {y, true, false}, communicator());
}

std::map<size_t, DistrArray::value_type> DistrArray::select_max_dot(size_t n, const DistrArray::SparseArray& y) const {
  auto name = std::string("DistrArray::select_max_dot:");
  if (y.empty())
    return {};
  if (size() < y.rbegin()->first + 1)
    error(name + " sparse array x is too large");
  if (n > size() || n > y.size())
    error(" n is too large");
  auto shifted_local_selection = std::map<size_t, value_type>();
  for_each_local_chunk({{*this, true, false}},
                       [n, &y, &shifted_local_selection](index_type start, size_t m, value_type* const* p) {
                         for (auto it = y.lower_bound(start); it != y.end() and it->first < start + m; ++it)
                           shifted_local_selection.emplace(it->first, std::abs(p[0][it->first - start] * it->second));
                         keep_largest(shifted_local_selection, n);
                       });
  return util::select_max_dot_broadcast(n, shifted_local_selection, communicator());
}

std::map<size_t, DistrArray::value_type> DistrArray::select(size_t n, bool max, bool ignore_sign) const {
  if (n > size())
    error("DistrArray::select: n is too large");
  // the selection is kept with values ordered so that larger is better
  auto shifted_local_selection = std::map<size_t, value_type>();
  for_each_local_chunk({{*this, true, false}},
                       [n, max, ignore_sign, &shifted_local_selection](index_type start, size_t m,
                                                                       value_type* const* p) {
                         const auto xs = Span<value_type>(p[0], m);
                         for (const auto& el : util::select<Span<value_type>, value_type>(std::min(n, m), xs, max,
                                                                                          ignore_sign))
                           shifted_local_selection.emplace(start + el.first, max ? el.second : -el.second);
                         keep_largest(shifted_local_selection, n);
                       });
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
  const auto [lo, hi] = local_range(x);
  const size_t length = hi - lo;
  auto nmin = length > size_t(n) ? size_t(n) : length;
  auto loc_extrema = std::list<std::pair<DistrArray::index_type, double>>();
  auto compare = Compare();
  auto compare_pair = [&compare](const auto& p1, const auto& p2) { return compare(p1.second, p2.second); };
  for_each_local_chunk({{x, true, false}}, [&, lo = lo](index_type start, size_t m, value_type* const* p) {
    for (size_t j = 0; j < m; ++j) {
      loc_extrema.emplace_back(start + j, p[0][j]);
      if (start + j - lo >= nmin) {
        loc_extrema.sort(compare_pair);
        loc_extrema.pop_back();
      }
    }
  });
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
    MPI_Igather(MPI_IN_PLACE, n, mpi_index_type(), indices_loc.data(), n, mpi_index_type(), 0, x.communicator(),
                &requests[1]);
    MPI_Igather(MPI_IN_PLACE, n, MPI_DOUBLE, values_loc.data(), n, MPI_DOUBLE, 0, x.communicator(),
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
    MPI_Igather(indices_loc.data(), n, mpi_index_type(), nullptr, n, mpi_index_type(), 0, x.communicator(),
                &requests[1]);
    MPI_Igather(values_loc.data(), n, MPI_DOUBLE, nullptr, n, MPI_DOUBLE, 0, x.communicator(), &requests[2]);
    MPI_Waitall(3, requests, MPI_STATUSES_IGNORE);
  }
  MPI_Ibcast(indices_glob.data(), n, mpi_index_type(), 0, x.communicator(), &requests[0]);
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
  array::copy({*this, false, true}, {y, true, false});
}

void DistrArray::copy_patch(const DistrArray& y, DistrArray::index_type start, DistrArray::index_type end) {
  auto name = std::string{"Array::copy_patch"};
  if (!compatible(y))
    error(name + " incompatible arrays");
  if (start > end)
    return;
  // offsets s to e within the local section, as calculated before the operation was paged
  const auto [lo, hi] = local_range(*this);
  const size_t length = hi - lo;
  const size_t s = start <= lo ? 0 : start - lo;
  const size_t e = end - start + 1 >= length ? length : end - start + 1;
  // the target is read as well as written, since only part of each page is copied
  for_each_local_chunk({{*this, true, true}, {y, true, false}},
                       [s, e, lo = lo](index_type chunk_start, size_t n, value_type* const* p) {
                         const size_t offset = chunk_start - lo;
                         for (size_t i = std::max(s, offset); i < std::min(e, offset + n); ++i)
                           p[0][i - offset] = p[1][i - offset];
                       });
}

DistrArray::value_type DistrArray::dot(const SparseArray& y) const {
  auto name = std::string{"Array::dot SparseArray "};
  if (y.empty())
    return 0;
  if (size() < y.rbegin()->first + 1)
    error(name + " sparse array x is incompatible");
  double res = 0;
  for_each_local_chunk({{*this, true, false}}, [&y, &res](index_type start, size_t n, value_type* const* p) {
    for (auto it = y.lower_bound(start); it != y.end() and it->first < start + n; ++it)
      res += p[0][it->first - start] * it->second;
  });
#ifdef HAVE_MPI_H
  MPI_Allreduce(MPI_IN_PLACE, &res, 1, MPI_DOUBLE, MPI_SUM, communicator());
#endif
  return res;
}

void DistrArray::axpy(value_type a, const SparseArray& y) {
  auto name = std::string{"Array::axpy SparseArray"};
  if (a == 0 || y.empty())
    return;
  if (size() < y.rbegin()->first + 1)
    error(name + " sparse array x is incompatible");
  for_each_local_chunk({{*this, true, true}}, [a, &y](index_type start, size_t n, value_type* const* p) {
    for (auto it = y.lower_bound(start); it != y.end() and it->first < start + n; ++it)
      p[0][it->first - start] += a * it->second;
  });
}

void DistrArray::copy(const SparseArray& y) {
  if (not y.empty() and size() < y.rbegin()->first + 1)
    error("Array::copy SparseArray sparse array y is incompatible");
  // write-only, so that an array on disk is written once and not read
  for_each_local_chunk({{*this, false, true}}, [&y](index_type start, size_t n, value_type* const* p) {
    std::fill_n(p[0], n, 0);
    for (auto it = y.lower_bound(start); it != y.end() and it->first < start + n; ++it)
      p[0][it->first - start] = it->second;
  });
}

void DistrArray::check_compatible(Span<const value_type> y, const std::string& name) const {
  if (y.size() != size())
    error(name + ": contiguous array of size " + std::to_string(y.size()) +
          " is incompatible with distributed array of size " + std::to_string(size()));
}

void DistrArray::copy(Span<const value_type> y) {
  check_compatible(y, "DistrArray::copy");
  array::copy({*this, false, true}, {y.data(), true, false});
}

void DistrArray::copy_to(Span<value_type> y) const {
  check_compatible(y, "DistrArray::copy_to");
  array::copy({y.data(), false, true}, {*this, true, false});
  util::gather_all(distribution(), communicator(), y.data());
}

void DistrArray::axpy(value_type a, Span<const value_type> y) {
  check_compatible(y, "DistrArray::axpy");
  array::axpy(a, {*this, true, true}, {y.data(), true, false});
}

void DistrArray::axpy_to(value_type a, Span<value_type> y) const {
  check_compatible(y, "DistrArray::axpy_to");
  if (a == 0)
    return;
  array::axpy(a, {y.data(), true, true}, {*this, true, false});
  util::gather_all(distribution(), communicator(), y.data());
}

DistrArray::value_type DistrArray::dot(Span<const value_type> y) const {
  auto prof = molpro::Profiler::single()->push("DistrArray::dot(Span)");
  check_compatible(y, "DistrArray::dot");
  return array::dot({*this, true, false}, {y.data(), true, false}, communicator());
}

std::map<size_t, DistrArray::value_type> DistrArray::select_max_dot(size_t n, Span<const value_type> y) const {
  check_compatible(y, "DistrArray::select_max_dot");
  if (n > size())
    error("DistrArray::select_max_dot: n is too large");
  return array::select_max_dot(n, {*this, true, false}, {y.data(), true, false}, communicator());
}

} // namespace molpro::linalg::array