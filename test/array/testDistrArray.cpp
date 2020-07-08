#include <algorithm>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <chrono>
#include <ga-mpi.h>
#include <ga.h>
#include <numeric>
#include <thread>

#include "parallel_util.h"
#include <molpro/linalg/array/DistrArrayGA.h>
#include <molpro/linalg/array/DistrArrayMPI3.h>

using ::testing::ContainerEq;
using ::testing::DoubleEq;
using ::testing::Each;
using ::testing::Pointwise;

using molpro::linalg::array::DistrArrayMPI3;
using molpro::linalg::array::util::LockMPI3;
using molpro::linalg::test::mpi_comm;

namespace {

int comm_rank(MPI_Comm comm) {
  int res;
  MPI_Comm_rank(comm, &res);
  return res;
}
int comm_size(MPI_Comm comm) {
  int res;
  MPI_Comm_size(comm, &res);
  return res;
}
} // namespace

TEST(LockMPI3, creation) { LockMPI3 lock{mpi_comm}; }

TEST(LockMPI3, lock_unlock) {
  LockMPI3 lock{mpi_comm};
  ASSERT_NO_FATAL_FAILURE(lock.lock());
  ASSERT_NO_FATAL_FAILURE(lock.unlock());
}

TEST(LockMPI3, scope) {
  LockMPI3 lock{mpi_comm};
  auto proxy = lock.scope();
}

TEST(LockMPI3, scope_and_delete) {
  auto lock = std::make_shared<LockMPI3>(mpi_comm);
  auto proxy = lock->scope();
  lock.reset();
}

// Simulate Fetch_and_op with put and get. If the lock mechanism does not work this is very likely to fail
// although it is not guranteed
TEST(LockMPI3, locking_mechanism_simulates_fetch_and_op) {
  auto comm = mpi_comm;
  auto size = comm_size(comm);
  auto rank = comm_rank(comm);
  if (size == 1)
    return;
  LockMPI3 lock{comm};
  MPI_Win win = MPI_WIN_NULL;
  int *base = nullptr;
  MPI_Win_allocate(1, sizeof(int), MPI_INFO_NULL, comm, &base, &win);
  MPI_Win_fence(0, win);
  int zero = 0;
  MPI_Put(&zero, 1, MPI_INT, rank, 0, 1, MPI_INT, win);
  MPI_Win_fence(0, win);
  MPI_Barrier(comm);
  int v = -1;
  {
    auto proxy = lock.scope();
    MPI_Get(&v, 1, MPI_INT, 0, 0, 1, MPI_INT, win);
    ++v;
    MPI_Put(&v, 1, MPI_INT, 0, 0, 1, MPI_INT, win);
  }
  MPI_Win_fence(0, win);
  if (rank == 0)
    v = base[0];
  MPI_Bcast(&v, 1, MPI_INT, 0, comm);
  {
    auto proxy = lock.scope();
    ASSERT_EQ(v, size);
  }
  MPI_Win_free(&win);
}

TEST(Array, constructor) {
  LockMPI3 lock{mpi_comm};
  auto proxy = lock.scope();
  size_t dim = 100;
  DistrArrayMPI3 a{dim, mpi_comm};
}

TEST(Array, constructor_copy) {
  LockMPI3 lock{mpi_comm};
  auto proxy = lock.scope();
  size_t dim = 100;
  DistrArrayMPI3 a{dim, mpi_comm};
  DistrArrayMPI3 b(a);
}

class ArrayInitializationF : public ::testing::Test, public DistrArrayMPI3 {
public:
  ArrayInitializationF()
      : DistrArrayMPI3(dim, mpi_comm), lock(m_communicator), m_comm_rank{comm_rank(m_communicator)} {};
  LockMPI3 lock;
  const int m_comm_rank;
  static const size_t dim = 30;
};

TEST_F(ArrayInitializationF, size) {
  {
    auto l = lock.scope();
    ASSERT_EQ(m_dimension, size());
  }
  sync();
}

TEST_F(ArrayInitializationF, empty) {
  {
    auto l = lock.scope();
    EXPECT_FALSE(m_allocated);
    ASSERT_TRUE(empty());
  }
  sync();
}

TEST_F(ArrayInitializationF, allocate_buffer) {
  allocate_buffer();
  {
    auto l = lock.scope();
    EXPECT_TRUE(m_allocated);
    ASSERT_FALSE(empty());
  }
  sync();
}

TEST_F(ArrayInitializationF, zero) {
  allocate_buffer();
  zero();
  sync();
  {
    auto l = lock.scope();
    EXPECT_TRUE(m_allocated);
    ASSERT_FALSE(empty());
  }
  sync();
}

TEST_F(ArrayInitializationF, vec) {
  allocate_buffer();
  zero();
  sync();
  {
    auto l = lock.scope();
    auto data = vec();
    auto empty_vec = std::vector<double>(data.size(), 0.);
    ASSERT_THAT(empty_vec, Pointwise(DoubleEq(), empty_vec));
  }
  sync();
}

TEST_F(ArrayInitializationF, get) {
  allocate_buffer();
  zero();
  sync();
  {
    auto l = lock.scope();
    auto data = get(0, 10);
    auto ref_vec = std::vector<double>(data.size(), 0.);
    ASSERT_THAT(ref_vec, Pointwise(DoubleEq(), data));
    data.assign(data.size(), 0);
    get(0, 10, data.data());
    ASSERT_THAT(ref_vec, Pointwise(DoubleEq(), data));
    data = get(0, dim - 1);
    auto same_as_vec = vec();
    ASSERT_THAT(data, Pointwise(DoubleEq(), same_as_vec));
  }
  sync();
}

TEST_F(ArrayInitializationF, put) {
  allocate_buffer();
  {
    auto l = lock.scope();
    auto range = std::vector<double>(dim);
    std::iota(range.begin(), range.end(), m_comm_rank);
    put(0, dim - 1, range.data());
    auto from_ga_buffer = get(0, dim - 1);
    ASSERT_THAT(from_ga_buffer, Pointwise(DoubleEq(), range));
  }
  sync();
}

TEST_F(ArrayInitializationF, fill) {
  allocate_buffer();
  zero();
  sync();
  auto ref_values = std::vector<double>(dim, 0.);
  auto from_ga_buffer = vec();
  {
    auto l = lock.scope();
    ASSERT_THAT(from_ga_buffer, Pointwise(DoubleEq(), ref_values));
  }
  sync();
  fill(42.0);
  sync();
  from_ga_buffer = vec();
  {
    auto l = lock.scope();
    std::fill(ref_values.begin(), ref_values.end(), 42.0);
    ASSERT_THAT(from_ga_buffer, Pointwise(DoubleEq(), ref_values));
  }
  sync();
}

class ArrayRangeF : public ::testing::Test, public DistrArrayMPI3 {
public:
  //! Stores a range in the buffer {1, 2, 3, 4, 5, .., dim}
  ArrayRangeF()
      : DistrArrayMPI3((size_t)dim, mpi_comm), lock(m_communicator), p_rank(comm_rank(m_communicator)),
        p_size(comm_size(m_communicator)) {
    DistrArrayMPI3::allocate_buffer();
    values.resize(dim);
    std::iota(values.begin(), values.end(), 1.);
    auto buffer = DistrArrayMPI3::local_buffer();
    std::copy(values.begin() + buffer->lo, values.begin() + buffer->hi, buffer->begin());
    sub_indices = {0, 1, 7, 15, 21, 29};
    for (auto el : sub_indices)
      sub_values.push_back(values[el]);
    DistrArrayMPI3::sync();
  }

  LockMPI3 lock;
  static const int dim = 30;
  int p_rank, p_size;
  std::vector<value_type> values;
  std::vector<index_type> sub_indices;
  std::vector<value_type> sub_values;
};

TEST_F(ArrayRangeF, gather) {
  sync();
  auto from_ga_buffer = gather(sub_indices);
  {
    auto l = lock.scope();
    ASSERT_THAT(from_ga_buffer, Pointwise(DoubleEq(), sub_values));
  }
  sync();
}

TEST_F(ArrayRangeF, scatter) {
  zero();
  sync();
  {
    auto proxy = lock.scope();
    scatter(sub_indices, sub_values);
    auto from_ga_buffer = gather(sub_indices);
    ASSERT_THAT(from_ga_buffer, Pointwise(DoubleEq(), sub_values));
    auto zero_values = std::vector<double>(sub_values.size(), 0.);
    scatter(sub_indices, zero_values);
  }
  sync();
}

TEST_F(ArrayRangeF, scatter_acc) {
  sync();
  {
    auto proxy = lock.scope();
    scatter_acc(sub_indices, sub_values);
    auto from_ga_buffer = gather(sub_indices);
    auto ref_values = sub_values;
    for (auto &el : ref_values)
      el *= 2;
    ASSERT_THAT(from_ga_buffer, Pointwise(DoubleEq(), ref_values));
    for (auto &el : sub_values)
      el *= -1;
    scatter_acc(sub_indices, sub_values);
  }
  sync();
}

TEST_F(ArrayRangeF, at) {
  sync();
  auto from_ga_buffer = std::vector<double>();
  for (int i = 0; i < dim; ++i) {
    from_ga_buffer.push_back(this->at(i));
  }
  {
    auto proxy = lock.scope();
    ASSERT_THAT(from_ga_buffer, Pointwise(DoubleEq(), values));
  }
  sync();
}

TEST_F(ArrayRangeF, minlocN) {
  int n = 10;
  auto ref_minloc_ind = std::vector<size_t>(n);
  std::iota(ref_minloc_ind.begin(), ref_minloc_ind.end(), 0);
  auto minloc_ind = min_loc_n(n);
  {
    auto proxy = lock.scope();
    ASSERT_THAT(minloc_ind, ContainerEq(ref_minloc_ind));
  }
  sync();
}

TEST_F(ArrayRangeF, minlocN_reverse) {
  int n = 10;
  std::reverse(values.begin(), values.end());
  if (p_rank == 0)
    put(0, dim - 1, values.data());
  sync();
  auto ref_minloc_ind = std::vector<index_type>(n);
  std::iota(ref_minloc_ind.rbegin(), ref_minloc_ind.rend(), dim - n);
  auto minloc_ind = min_loc_n(n);
  {
    auto proxy = lock.scope();
    ASSERT_THAT(minloc_ind, ContainerEq(ref_minloc_ind));
  }
  sync();
}

TEST_F(ArrayRangeF, max_n) {
  int n = 10;
  auto ref_minloc_ind = std::vector<size_t>(n);
  std::iota(ref_minloc_ind.rbegin(), ref_minloc_ind.rend(), size() - n);
  auto maxloc = max_n(n);
  auto maxloc_ind = std::vector<size_t>(n);
  std::transform(maxloc.cbegin(), maxloc.cend(), maxloc_ind.begin(), [](const auto &p) { return p.first; });
  {
    auto proxy = lock.scope();
    ASSERT_THAT(maxloc_ind, ContainerEq(ref_minloc_ind));
  }
  sync();
}

TEST_F(ArrayRangeF, min_abs_n) {
  auto loc_buffer = local_buffer();
  for (size_t i = 0; i < loc_buffer->size(); ++i)
    if ((i + loc_buffer->lo) % 2 == 1)
      loc_buffer->at(i) *= -1;
  int n = 10;
  auto ref_minloc_ind = std::vector<size_t>(n);
  std::iota(ref_minloc_ind.begin(), ref_minloc_ind.end(), 0);
  sync();
  auto minloc = min_abs_n(n);
  auto minloc_ind = std::vector<size_t>(n);
  std::transform(minloc.cbegin(), minloc.cend(), minloc_ind.begin(), [](const auto &p) { return p.first; });
  {
    auto proxy = lock.scope();
    ASSERT_THAT(minloc_ind, ContainerEq(ref_minloc_ind));
  }
  sync();
}

TEST_F(ArrayRangeF, max_abs_n) {
  auto loc_buffer = local_buffer();
  for (size_t i = 0; i < loc_buffer->size(); ++i)
    if ((i + loc_buffer->lo) % 2 == 1)
      loc_buffer->at(i) *= -1;
  int n = 10;
  auto ref_minloc_ind = std::vector<size_t>(n);
  std::iota(ref_minloc_ind.rbegin(), ref_minloc_ind.rend(), size() - n);
  sync();
  auto maxloc = max_abs_n(n);
  auto maxloc_ind = std::vector<size_t>(n);
  std::transform(maxloc.cbegin(), maxloc.cend(), maxloc_ind.begin(), [](const auto &p) { return p.first; });
  {
    auto proxy = lock.scope();
    ASSERT_THAT(maxloc_ind, ContainerEq(ref_minloc_ind));
  }
  sync();
}

TEST_F(ArrayRangeF, scal_double) {
  double alpha = 1.5;
  scal(alpha);
  sync();
  for (auto &el : values)
    el *= alpha;
  auto from_ga_buffer = vec();
  {
    auto proxy = lock.scope();
    ASSERT_THAT(from_ga_buffer, Pointwise(DoubleEq(), values));
  }
  sync();
}

TEST_F(ArrayRangeF, add_double) {
  double alpha = 100.;
  add(alpha);
  sync();
  for (auto &el : values)
    el += alpha;
  {
    auto proxy = lock.scope();
    auto from_ga_buffer = vec();
    ASSERT_THAT(from_ga_buffer, Pointwise(DoubleEq(), values));
  }
  sync();
}

TEST_F(ArrayRangeF, sub_double) {
  double alpha = 100.;
  sub(alpha);
  sync();
  for (auto &el : values)
    el -= alpha;
  {
    auto proxy = lock.scope();
    auto from_ga_buffer = vec();
    ASSERT_THAT(from_ga_buffer, Pointwise(DoubleEq(), values));
  }
  sync();
}

TEST_F(ArrayRangeF, recip) {
  recip();
  sync();
  for (auto &el : values)
    el = 1. / el;
  {
    auto proxy = lock.scope();
    auto from_ga_buffer = vec();
    ASSERT_THAT(from_ga_buffer, Pointwise(DoubleEq(), values));
  }
  sync();
}

class ArrayCollectiveOpF : public ::testing::Test {
public:
  ArrayCollectiveOpF() : lock(mpi_comm), p_rank(GA_Nodeid()), p_size(GA_Nnodes()), a(dim, mpi_comm), b(dim, mpi_comm) {
    a.allocate_buffer();
    b.allocate_buffer();
    range_alpha.resize(dim);
    range_beta.resize(dim);
    std::iota(range_alpha.begin(), range_alpha.end(), 0);
    for (const auto &el : range_alpha) {
      range_beta.push_back(el * beta);
    }
    auto n = dim / every_nth_element;
    auto indices = std::vector<size_t>(n);
    for (int i = 0; i < n; ++i) {
      indices[i] = range_beta[i * every_nth_element];
    }
    for (const auto &el : indices) {
      sparse_array.emplace(el, range_beta[el]);
    }
  }

  LockMPI3 lock;
  static const int dim = 30;
  const double alpha = 1;
  const double beta = 2;
  std::vector<double> range_alpha;
  std::vector<double> range_beta;
  int p_rank, p_size;
  DistrArrayMPI3 a;
  DistrArrayMPI3 b;
  static const int every_nth_element = 5; // period for selection of sparse array elements
  std::map<size_t, double> sparse_array;  // sparse selection of range_beta elements
};

TEST_F(ArrayCollectiveOpF, add) {
  a.fill(alpha);
  b.fill(beta);
  a.add(b);
  a.sync();
  auto from_ga_buffer_a = a.vec();
  auto from_ga_buffer_b = b.vec();
  {
    auto proxy = lock.scope();
    ASSERT_THAT(from_ga_buffer_a, Each(DoubleEq(alpha + beta)));
    ASSERT_THAT(from_ga_buffer_b, Each(DoubleEq(beta)));
  }
  a.sync();
}

TEST_F(ArrayCollectiveOpF, sub) {
  a.fill(alpha);
  b.fill(beta);
  a.sub(b);
  a.sync();
  auto from_ga_buffer_a = a.vec();
  auto from_ga_buffer_b = b.vec();
  {
    auto proxy = lock.scope();
    ASSERT_THAT(from_ga_buffer_a, Each(DoubleEq(alpha - beta)));
    ASSERT_THAT(from_ga_buffer_b, Each(DoubleEq(beta)));
  }
  a.sync();
}

TEST_F(ArrayCollectiveOpF, axpy) {
  double scale = -3.0;
  a.fill(alpha);
  b.fill(beta);
  a.axpy(scale, b);
  a.sync();
  auto from_ga_buffer_a = a.vec();
  auto from_ga_buffer_b = b.vec();
  {
    auto proxy = lock.scope();
    ASSERT_THAT(from_ga_buffer_a, Each(DoubleEq(alpha + scale * beta)));
    ASSERT_THAT(from_ga_buffer_b, Each(DoubleEq(beta)));
  }
  a.sync();
}

TEST_F(ArrayCollectiveOpF, axpy_map) {
  a.fill(alpha);
  double scale = 5.0;
  a.axpy(scale, sparse_array);
  a.sync();
  auto from_ga_buffer_a = a.vec();
  auto ref_vals = std::vector<double>(dim, alpha);
  for (const auto &item : sparse_array) {
    ref_vals[item.first] += scale * item.second;
  }
  {
    auto proxy = lock.scope();
    ASSERT_THAT(from_ga_buffer_a, Pointwise(DoubleEq(), ref_vals));
  }
  a.sync();
}

TEST_F(ArrayCollectiveOpF, dot_Array) {
  if (p_rank == 0) {
    a.put(0, dim - 1, range_alpha.data());
    b.put(0, dim - 1, range_beta.data());
  }
  a.sync();
  auto ga_dot = a.dot(b);
  double ref_dot = std::inner_product(range_alpha.begin(), range_alpha.end(), range_beta.begin(), 0.);
  {
    auto proxy = lock.scope();
    ASSERT_THAT(ga_dot, DoubleEq(ref_dot));
  }
  a.sync();
}

TEST_F(ArrayCollectiveOpF, dot_map) {
  if (p_rank == 0)
    a.put(0, dim - 1, range_alpha.data());
  a.sync();
  auto ga_dot = a.dot(sparse_array);
  double ref_dot = 0.;
  for (const auto &item : sparse_array) {
    ref_dot += range_alpha[item.first] * item.second;
  }
  {
    auto proxy = lock.scope();
    ASSERT_THAT(ga_dot, DoubleEq(ref_dot));
  }
  a.sync();
}

TEST_F(ArrayCollectiveOpF, times) {
  a.fill(alpha);
  b.fill(beta);
  DistrArrayMPI3 c{a};
  c.times(a, b);
  a.sync();
  auto from_ga_buffer_a = a.vec();
  auto from_ga_buffer_b = b.vec();
  auto from_ga_buffer_c = c.vec();
  auto ref_a = std::vector<double>(dim, alpha);
  auto ref_b = std::vector<double>(dim, beta);
  auto ref_c = std::vector<double>(dim, alpha * beta);
  {
    auto proxy = lock.scope();
    ASSERT_THAT(from_ga_buffer_a, Pointwise(DoubleEq(), ref_a));
    ASSERT_THAT(from_ga_buffer_b, Pointwise(DoubleEq(), ref_b));
    ASSERT_THAT(from_ga_buffer_c, Pointwise(DoubleEq(), ref_c));
  }
  a.sync();
}

// c[i] -= a[i]/(b[i]+shift)
TEST_F(ArrayCollectiveOpF, divide_append_negative) {
  double shift = 0.5;
  a.fill(alpha);
  b.fill(beta);
  DistrArrayMPI3 c{a};
  c.divide(a, b, shift, true, true);
  a.sync();
  auto from_ga_buffer_a = a.vec();
  auto from_ga_buffer_b = b.vec();
  auto from_ga_buffer_c = c.vec();
  auto ref_a = std::vector<double>(dim, alpha);
  auto ref_b = std::vector<double>(dim, beta);
  auto ref_c = std::vector<double>(dim, alpha - alpha / (beta + shift));
  {
    auto proxy = lock.scope();
    ASSERT_THAT(from_ga_buffer_a, Pointwise(DoubleEq(), ref_a));
    ASSERT_THAT(from_ga_buffer_b, Pointwise(DoubleEq(), ref_b));
    ASSERT_THAT(from_ga_buffer_c, Pointwise(DoubleEq(), ref_c));
  }
  a.sync();
}

// c[i] += a[i]/(b[i]+shift)
TEST_F(ArrayCollectiveOpF, divide_append_positive) {
  double shift = 0.5;
  a.fill(alpha);
  b.fill(beta);
  DistrArrayMPI3 c{a};
  c.divide(a, b, shift, true, false);
  a.sync();
  auto from_ga_buffer_a = a.vec();
  auto from_ga_buffer_b = b.vec();
  auto from_ga_buffer_c = c.vec();
  auto ref_a = std::vector<double>(dim, alpha);
  auto ref_b = std::vector<double>(dim, beta);
  auto ref_c = std::vector<double>(dim, alpha + alpha / (beta + shift));
  {
    auto proxy = lock.scope();
    ASSERT_THAT(from_ga_buffer_a, Pointwise(DoubleEq(), ref_a));
    ASSERT_THAT(from_ga_buffer_b, Pointwise(DoubleEq(), ref_b));
    ASSERT_THAT(from_ga_buffer_c, Pointwise(DoubleEq(), ref_c));
  }
  a.sync();
}

// c[i] = a[i]/(b[i]+shift)
TEST_F(ArrayCollectiveOpF, divide_overwrite_positive) {
  double shift = 0.5;
  a.fill(alpha);
  b.fill(beta);
  DistrArrayMPI3 c{a};
  c.divide(a, b, shift, false, false);
  a.sync();
  auto from_ga_buffer_a = a.vec();
  auto from_ga_buffer_b = b.vec();
  auto from_ga_buffer_c = c.vec();
  auto ref_a = std::vector<double>(dim, alpha);
  auto ref_b = std::vector<double>(dim, beta);
  auto ref_c = std::vector<double>(dim, alpha / (beta + shift));
  {
    auto proxy = lock.scope();
    ASSERT_THAT(from_ga_buffer_a, Pointwise(DoubleEq(), ref_a));
    ASSERT_THAT(from_ga_buffer_b, Pointwise(DoubleEq(), ref_b));
    ASSERT_THAT(from_ga_buffer_c, Pointwise(DoubleEq(), ref_c));
  }
  a.sync();
}

struct DistributionF : public ::testing::Test, public DistrArrayMPI3 {
  DistributionF() : DistrArrayMPI3{m_dim, mpi_comm} {}

  DistributionMPI3 factory(int n_proc, size_t dim) { return DistributionMPI3(n_proc, dim); }
  static const size_t m_dim = 30;
};

TEST_F(DistributionF, constructor) {
  int np = 5;
  auto d = factory(np, 17);
  auto ref_distribution = std::vector<std::pair<unsigned long, size_t>>{{0, 4}, {4, 4}, {8, 3}, {11, 3}, {14, 3}};
  std::vector<decltype((d.range)(1))> proc_buffer;
  for (size_t i = 0; i < np; ++i)
    proc_buffer.emplace_back(d.range(i));
  ASSERT_THAT(proc_buffer, ContainerEq(ref_distribution));
}

TEST_F(DistributionF, process_map) {
  int np = 3;
  auto d = factory(np, 20);
  auto ref_lo_hi = std::vector<std::pair<unsigned long, unsigned long>>{{0, 0},   {1, 5},   {3, 7},   {5, 9},
                                                                        {13, 14}, {15, 19}, {18, 12}, {17, 6}};
  auto ref_process_pairs =
      std::vector<std::pair<int, int>>{{0, 0}, {0, 0}, {0, 1}, {0, 1}, {1, 2}, {2, 2}, {2, 1}, {2, 0}};
  auto process_pairs = std::vector<std::pair<int, int>>(ref_lo_hi.size());
  std::transform(ref_lo_hi.begin(), ref_lo_hi.end(), process_pairs.begin(),
                 [&d](auto p) { return d.locate_process(p.first, p.second); });
  ASSERT_THAT(process_pairs, ContainerEq(ref_process_pairs));
}
