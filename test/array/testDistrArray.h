#ifndef LINEARALGEBRA_TEST_ARRAY_TESTDISTRARRAY_H
#define LINEARALGEBRA_TEST_ARRAY_TESTDISTRARRAY_H
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <functional>
#include <numeric>
#include <thread>

#include "parallel_util.h"
#include <molpro/linalg/array/util.h>

using ::testing::ContainerEq;
using ::testing::DoubleEq;
using ::testing::Each;
using ::testing::Pointwise;

using molpro::linalg::array::util::LockMPI3;
using molpro::linalg::test::mpi_comm;

template <class Array>
struct TestDistrArray : public ::testing::Test {};

TYPED_TEST_SUITE_P(TestDistrArray);

TYPED_TEST_P(TestDistrArray, constructor) {
  size_t dim = 100;
  TypeParam a{dim, mpi_comm};
}

TYPED_TEST_P(TestDistrArray, constructor_copy) {
  LockMPI3 lock{mpi_comm};
  size_t dim = 100;
  TypeParam a{dim, mpi_comm};
  TypeParam b(a);
  auto proxy = lock.scope();
  ASSERT_EQ(b.size(), a.size());
  ASSERT_EQ(b.communicator(), a.communicator());
}

TYPED_TEST_P(TestDistrArray, constructor_copy_allocated) {
  LockMPI3 lock{mpi_comm};
  size_t dim = 100;
  double alpha = 1;
  TypeParam a{dim, mpi_comm};
  a.fill(alpha);
  TypeParam b(a);
  auto proxy = lock.scope();
  ASSERT_EQ(b.size(), a.size());
  ASSERT_EQ(b.communicator(), a.communicator());
  ASSERT_THAT(*b.local_buffer(), Each(DoubleEq(alpha)));
}

TYPED_TEST_P(TestDistrArray, copy_assignment_op_allocated) {
  LockMPI3 lock{mpi_comm};
  size_t dim = 100;
  double alpha = 1;
  double beta = 11;
  TypeParam a{dim, mpi_comm};
  a.fill(alpha);
  TypeParam b{dim, mpi_comm};
  b = a;
  auto proxy = lock.scope();
  ASSERT_EQ(b.size(), a.size());
  ASSERT_EQ(b.communicator(), a.communicator());
  ASSERT_THAT(*b.local_buffer(), Each(DoubleEq(alpha)));
}

TYPED_TEST_P(TestDistrArray, constructor_move) {
  LockMPI3 lock{mpi_comm};
  size_t dim = 100;
  TypeParam &&a{dim, mpi_comm};
  TypeParam b(std::move(a));
  auto proxy = lock.scope();
  ASSERT_EQ(b.size(), dim);
  ASSERT_EQ(b.communicator(), mpi_comm);
}

TYPED_TEST_P(TestDistrArray, constructor_move_allocated) {
  LockMPI3 lock{mpi_comm};
  size_t dim = 100;
  double alpha = 1;
  TypeParam &&a{dim, mpi_comm};
  a.fill(alpha);
  TypeParam b(std::move(a));
  auto proxy = lock.scope();
  ASSERT_EQ(b.size(), dim);
  ASSERT_EQ(b.communicator(), mpi_comm);
  ASSERT_THAT(*b.local_buffer(), Each(DoubleEq(alpha)));
}

TYPED_TEST_P(TestDistrArray, move_assignment_op_allocated) {
  LockMPI3 lock{mpi_comm};
  size_t dim = 100;
  double alpha = 1;
  double beta = 11;
  TypeParam &&a{dim, mpi_comm};
  a.fill(alpha);
  TypeParam b{dim, mpi_comm};
  b.fill(beta);
  b = std::move(a);
  auto proxy = lock.scope();
  ASSERT_EQ(b.size(), dim);
  ASSERT_EQ(b.communicator(), mpi_comm);
  ASSERT_THAT(*b.local_buffer(), Each(DoubleEq(alpha)));
}

TYPED_TEST_P(TestDistrArray, select_max_dot) {
  LockMPI3 lock{mpi_comm};
  const size_t dim = 30;
  const size_t n = 5; // values to select
  TypeParam x{dim, mpi_comm};
  TypeParam y{dim, mpi_comm};
  y.fill(1);
  int rank, comm_size;
  MPI_Comm_rank(mpi_comm, &rank);
  MPI_Comm_size(mpi_comm, &comm_size);
  auto buffer = x.local_buffer();
  std::iota(buffer->begin(), buffer->end(), typename TypeParam::value_type(buffer->start()));
  x.sync();
  auto selection_ref = std::map<size_t, typename TypeParam::value_type>();
  for (size_t i = 0; i < n; ++i)
    selection_ref[dim - 1 - i] = dim - 1 - i;
  auto selection = x.select_max_dot(n, y);
  {
    auto l = lock.scope();
    ASSERT_THAT(selection, ContainerEq(selection_ref));
  }
}

template <typename Array>
class DistArrayInitialization : public Array {
public:
  DistArrayInitialization() : Array(dim, mpi_comm), lock(mpi_comm), m_comm_rank{0} {
    MPI_Comm_rank(mpi_comm, &m_comm_rank);
  };
  LockMPI3 lock;
  int m_comm_rank;
  static const size_t dim = 30;
};

template <typename Array>
class DistArrayBasicF : public DistArrayInitialization<Array>, public ::testing::Test {};

TYPED_TEST_SUITE_P(DistArrayBasicF);

TYPED_TEST_P(DistArrayBasicF, size) {
  {
    auto l = this->lock.scope();
    ASSERT_EQ(TypeParam::m_dimension, TypeParam::size());
  }
  TypeParam::sync();
}

TYPED_TEST_P(DistArrayBasicF, zero) {
  TypeParam::zero();
  TypeParam::sync();
}

template <typename Array>
class DistArrayBasicRMAF : public DistArrayInitialization<Array>, public ::testing::Test {};

TYPED_TEST_SUITE_P(DistArrayBasicRMAF);

TYPED_TEST_P(DistArrayBasicRMAF, vec) {
  TypeParam::zero();
  TypeParam::sync();
  {
    auto l = this->lock.scope();
    auto data = TypeParam::vec();
    auto empty_vec = std::vector<double>(data.size(), 0.);
    ASSERT_THAT(empty_vec, Pointwise(DoubleEq(), empty_vec));
  }
  TypeParam::sync();
}

TYPED_TEST_P(DistArrayBasicRMAF, get) {
  TypeParam::zero();
  TypeParam::sync();
  {
    auto l = this->lock.scope();
    auto data = TypeParam::get(0, 10);
    auto ref_vec = std::vector<double>(data.size(), 0.);
    ASSERT_THAT(ref_vec, Pointwise(DoubleEq(), data));
    data.assign(data.size(), 0);
    TypeParam::get(0, 10, data.data());
    ASSERT_THAT(ref_vec, Pointwise(DoubleEq(), data));
    data = TypeParam::get(0, this->dim);
    auto same_as_vec = TypeParam::vec();
    ASSERT_THAT(data, Pointwise(DoubleEq(), same_as_vec));
  }
  TypeParam::sync();
}

TYPED_TEST_P(DistArrayBasicRMAF, put) {
  {
    auto l = this->lock.scope();
    auto range = std::vector<double>(this->dim);
    std::iota(range.begin(), range.end(), this->m_comm_rank);
    TypeParam::put(0, this->dim, range.data());
    auto from_ga_buffer = TypeParam::get(0, this->dim);
    ASSERT_THAT(from_ga_buffer, Pointwise(DoubleEq(), range));
  }
  TypeParam::sync();
}

TYPED_TEST_P(DistArrayBasicF, fill) {
  TypeParam::zero();
  TypeParam::sync();
  auto ref_values = std::vector<double>(this->dim, 0.);
  auto from_ga_buffer = TypeParam::vec();
  {
    auto l = this->lock.scope();
    ASSERT_THAT(from_ga_buffer, Pointwise(DoubleEq(), ref_values));
  }
  TypeParam::sync();
  TypeParam::fill(42.0);
  TypeParam::sync();
  from_ga_buffer = TypeParam::vec();
  {
    auto l = this->lock.scope();
    std::fill(ref_values.begin(), ref_values.end(), 42.0);
    ASSERT_THAT(from_ga_buffer, Pointwise(DoubleEq(), ref_values));
  }
  TypeParam::sync();
}

template <typename Array>
class DistrArrayRangeF : public Array {
public:
  //! Stores a range in the buffer {1, 2, 3, 4, 5, .., dim}
  DistrArrayRangeF() : Array((size_t)dim, mpi_comm), lock(mpi_comm), p_rank(0), p_size(0) {
    MPI_Comm_rank(mpi_comm, &p_rank);
    MPI_Comm_size(mpi_comm, &p_size);
    values.resize(dim);
    std::iota(values.begin(), values.end(), 1.);
    auto buffer = Array::local_buffer();
    std::copy(values.begin() + buffer->start(), values.begin() + buffer->start() + buffer->size(), buffer->begin());
    sub_indices = {0, 1, 7, 15, 21, 29};
    for (auto el : sub_indices)
      sub_values.push_back(values[el]);
    Array::sync();
  }

  LockMPI3 lock;
  static const int dim = 30;
  int p_rank, p_size;
  std::vector<typename Array::value_type> values;
  std::vector<typename Array::index_type> sub_indices;
  std::vector<typename Array::value_type> sub_values;
};

template <typename Array>
class DistrArrayRangeRMAF : public DistrArrayRangeF<Array>, public ::testing::Test {};

TYPED_TEST_SUITE_P(DistrArrayRangeRMAF);

TYPED_TEST_P(DistrArrayRangeRMAF, gather) {
  TypeParam::sync();
  auto from_ga_buffer = TypeParam::gather(this->sub_indices);
  {
    auto l = this->lock.scope();
    ASSERT_THAT(from_ga_buffer, Pointwise(DoubleEq(), this->sub_values));
  }
  TypeParam::sync();
}

TYPED_TEST_P(DistrArrayRangeRMAF, scatter) {
  TypeParam::zero();
  TypeParam::sync();
  {
    auto proxy = this->lock.scope();
    TypeParam::scatter(this->sub_indices, this->sub_values);
    auto from_ga_buffer = TypeParam::gather(this->sub_indices);
    ASSERT_THAT(from_ga_buffer, Pointwise(DoubleEq(), this->sub_values));
    auto zero_values = std::vector<double>(this->sub_values.size(), 0.);
    TypeParam::scatter(this->sub_indices, zero_values);
  }
  TypeParam::sync();
}

TYPED_TEST_P(DistrArrayRangeRMAF, scatter_acc) {
  TypeParam::sync();
  {
    auto proxy = this->lock.scope();
    TypeParam::scatter_acc(this->sub_indices, this->sub_values);
    auto from_ga_buffer = TypeParam::gather(this->sub_indices);
    auto ref_values = this->sub_values;
    for (auto &el : ref_values)
      el *= 2;
    ASSERT_THAT(from_ga_buffer, Pointwise(DoubleEq(), ref_values));
    for (auto &el : this->sub_values)
      el *= -1;
    TypeParam::scatter_acc(this->sub_indices, this->sub_values);
  }
  TypeParam::sync();
}

TYPED_TEST_P(DistrArrayRangeRMAF, at) {
  TypeParam::sync();
  auto from_ga_buffer = std::vector<double>();
  for (int i = 0; i < this->dim; ++i) {
    from_ga_buffer.push_back(TypeParam::at(i));
  }
  {
    auto proxy = this->lock.scope();
    ASSERT_THAT(from_ga_buffer, Pointwise(DoubleEq(), this->values));
  }
  TypeParam::sync();
}

template <typename Array>
class DistrArrayRangeMinMaxF : public DistrArrayRangeF<Array>, public ::testing::Test {};

TYPED_TEST_SUITE_P(DistrArrayRangeMinMaxF);

TYPED_TEST_P(DistrArrayRangeMinMaxF, min_loc_n) {
  int n = 10;
  auto ref_minloc_ind = std::vector<size_t>(n);
  std::iota(ref_minloc_ind.begin(), ref_minloc_ind.end(), 0);
  auto minloc_ind = TypeParam::min_loc_n(n);
  {
    auto proxy = this->lock.scope();
    ASSERT_THAT(minloc_ind, ContainerEq(ref_minloc_ind));
  }
  TypeParam::sync();
}

TYPED_TEST_P(DistrArrayRangeMinMaxF, min_loc_n_reverse) {
  int n = 10;
  std::reverse(this->values.begin(), this->values.end());
  if (this->p_rank == 0)
    TypeParam::put(0, this->dim, this->values.data());
  TypeParam::sync();
  auto ref_minloc_ind = std::vector<typename TypeParam::index_type>(n);
  std::iota(ref_minloc_ind.rbegin(), ref_minloc_ind.rend(), this->dim - n);
  auto minloc_ind = TypeParam::min_loc_n(n);
  {
    auto proxy = this->lock.scope();
    ASSERT_THAT(minloc_ind, ContainerEq(ref_minloc_ind));
  }
  TypeParam::sync();
}

TYPED_TEST_P(DistrArrayRangeMinMaxF, max_n) {
  int n = 10;
  auto ref_minloc_ind = std::vector<size_t>(n);
  std::iota(ref_minloc_ind.rbegin(), ref_minloc_ind.rend(), TypeParam::size() - n);
  auto maxloc = TypeParam::max_n(n);
  auto maxloc_ind = std::vector<size_t>(n);
  std::transform(maxloc.cbegin(), maxloc.cend(), maxloc_ind.begin(), [](const auto &p) { return p.first; });
  {
    auto proxy = this->lock.scope();
    ASSERT_THAT(maxloc_ind, ContainerEq(ref_minloc_ind));
  }
  TypeParam::sync();
}

TYPED_TEST_P(DistrArrayRangeMinMaxF, min_abs_n) {
  auto loc_buffer = TypeParam::local_buffer();
  for (size_t i = 0; i < loc_buffer->size(); ++i)
    if ((i + loc_buffer->start()) % 2 == 1)
      (*loc_buffer)[i] *= -1;
  int n = 10;
  auto ref_minloc_ind = std::vector<size_t>(n);
  std::iota(ref_minloc_ind.begin(), ref_minloc_ind.end(), 0);
  TypeParam::sync();
  auto minloc = TypeParam::min_abs_n(n);
  auto minloc_ind = std::vector<size_t>(n);
  std::transform(minloc.cbegin(), minloc.cend(), minloc_ind.begin(), [](const auto &p) { return p.first; });
  {
    auto proxy = this->lock.scope();
    ASSERT_THAT(minloc_ind, ContainerEq(ref_minloc_ind));
  }
  TypeParam::sync();
}

TYPED_TEST_P(DistrArrayRangeMinMaxF, max_abs_n) {
  auto loc_buffer = TypeParam::local_buffer();
  for (size_t i = 0; i < loc_buffer->size(); ++i)
    if ((i + loc_buffer->start()) % 2 == 1)
      (*loc_buffer)[i] *= -1;
  int n = 10;
  auto ref_minloc_ind = std::vector<size_t>(n);
  std::iota(ref_minloc_ind.rbegin(), ref_minloc_ind.rend(), TypeParam::size() - n);
  TypeParam::sync();
  auto maxloc = TypeParam::max_abs_n(n);
  auto maxloc_ind = std::vector<size_t>(n);
  std::transform(maxloc.cbegin(), maxloc.cend(), maxloc_ind.begin(), [](const auto &p) { return p.first; });
  {
    auto proxy = this->lock.scope();
    ASSERT_THAT(maxloc_ind, ContainerEq(ref_minloc_ind));
  }
  TypeParam::sync();
}

template <typename Array>
class DistrArrayRangeLinAlgF : public DistrArrayRangeF<Array>, public ::testing::Test {};

TYPED_TEST_SUITE_P(DistrArrayRangeLinAlgF);

TYPED_TEST_P(DistrArrayRangeLinAlgF, scal_double) {
  double alpha = 1.5;
  TypeParam::scal(alpha);
  TypeParam::sync();
  for (auto &el : this->values)
    el *= alpha;
  auto from_ga_buffer = TypeParam::vec();
  {
    auto proxy = this->lock.scope();
    ASSERT_THAT(from_ga_buffer, Pointwise(DoubleEq(), this->values));
  }
  TypeParam::sync();
}

TYPED_TEST_P(DistrArrayRangeLinAlgF, add_double) {
  double alpha = 100.;
  TypeParam::add(alpha);
  TypeParam::sync();
  for (auto &el : this->values)
    el += alpha;
  {
    auto proxy = this->lock.scope();
    auto from_ga_buffer = TypeParam::vec();
    ASSERT_THAT(from_ga_buffer, Pointwise(DoubleEq(), this->values));
  }
  TypeParam::sync();
}

TYPED_TEST_P(DistrArrayRangeLinAlgF, sub_double) {
  double alpha = 100.;
  TypeParam::sub(alpha);
  sync();
  for (auto &el : this->values)
    el -= alpha;
  {
    auto proxy = this->lock.scope();
    auto from_ga_buffer = TypeParam::vec();
    ASSERT_THAT(from_ga_buffer, Pointwise(DoubleEq(), this->values));
  }
  TypeParam::sync();
}

TYPED_TEST_P(DistrArrayRangeLinAlgF, recip) {
  TypeParam::recip();
  TypeParam::sync();
  for (auto &el : this->values)
    el = 1. / el;
  {
    auto proxy = this->lock.scope();
    auto from_ga_buffer = TypeParam::vec();
    ASSERT_THAT(from_ga_buffer, Pointwise(DoubleEq(), this->values));
  }
  TypeParam::sync();
}

template <typename Array>
class DistrArrayCollectiveLinAlgF : public ::testing::Test {
public:
  DistrArrayCollectiveLinAlgF() : lock(mpi_comm), p_rank(0), p_size(0), a(dim, mpi_comm), b(dim, mpi_comm) {
    MPI_Comm_rank(mpi_comm, &p_rank);
    MPI_Comm_size(mpi_comm, &p_size);
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
  Array a;
  Array b;
  static const int every_nth_element = 5; // period for selection of sparse array elements
  std::map<size_t, double> sparse_array;  // sparse selection of range_beta elements
};

TYPED_TEST_SUITE_P(DistrArrayCollectiveLinAlgF);

TYPED_TEST_P(DistrArrayCollectiveLinAlgF, add) {
  this->a.fill(this->alpha);
  this->b.fill(this->beta);
  this->a.add(this->b);
  this->a.sync();
  auto from_ga_buffer_a = this->a.vec();
  auto from_ga_buffer_b = this->b.vec();
  {
    auto proxy = this->lock.scope();
    ASSERT_THAT(from_ga_buffer_a, Each(DoubleEq(this->alpha + this->beta)));
    ASSERT_THAT(from_ga_buffer_b, Each(DoubleEq(this->beta)));
  }
  this->a.sync();
}

TYPED_TEST_P(DistrArrayCollectiveLinAlgF, sub) {
  this->a.fill(this->alpha);
  this->b.fill(this->beta);
  this->a.sub(this->b);
  this->a.sync();
  auto from_ga_buffer_a = this->a.vec();
  auto from_ga_buffer_b = this->b.vec();
  {
    auto proxy = this->lock.scope();
    ASSERT_THAT(from_ga_buffer_a, Each(DoubleEq(this->alpha - this->beta)));
    ASSERT_THAT(from_ga_buffer_b, Each(DoubleEq(this->beta)));
  }
  this->a.sync();
}

TYPED_TEST_P(DistrArrayCollectiveLinAlgF, axpy) {
  double scale = -3.0;
  this->a.fill(this->alpha);
  this->b.fill(this->beta);
  this->a.axpy(scale, this->b);
  this->a.sync();
  auto from_ga_buffer_a = this->a.vec();
  auto from_ga_buffer_b = this->b.vec();
  {
    auto proxy = this->lock.scope();
    ASSERT_THAT(from_ga_buffer_a, Each(DoubleEq(this->alpha + scale * this->beta)));
    ASSERT_THAT(from_ga_buffer_b, Each(DoubleEq(this->beta)));
  }
  this->a.sync();
}

TYPED_TEST_P(DistrArrayCollectiveLinAlgF, axpy_map) {
  this->a.fill(this->alpha);
  double scale = 5.0;
  this->a.axpy(scale, this->sparse_array);
  this->a.sync();
  auto from_ga_buffer_a = this->a.vec();
  auto ref_vals = std::vector<double>(this->dim, this->alpha);
  for (const auto &item : this->sparse_array) {
    ref_vals[item.first] += scale * item.second;
  }
  {
    auto proxy = this->lock.scope();
    ASSERT_THAT(from_ga_buffer_a, Pointwise(DoubleEq(), ref_vals));
  }
  this->a.sync();
}

TYPED_TEST_P(DistrArrayCollectiveLinAlgF, dot_array) {
  if (this->p_rank == 0) {
    this->a.put(0, this->dim, this->range_alpha.data());
    this->b.put(0, this->dim, this->range_beta.data());
  }
  this->a.sync();
  auto ga_dot = this->a.dot(this->b);
  double ref_dot = std::inner_product(this->range_alpha.begin(), this->range_alpha.end(), this->range_beta.begin(), 0.);
  {
    auto proxy = this->lock.scope();
    ASSERT_THAT(ga_dot, DoubleEq(ref_dot));
  }
  this->a.sync();
}

TYPED_TEST_P(DistrArrayCollectiveLinAlgF, dot_map) {
  if (this->p_rank == 0)
    this->a.put(0, this->dim, this->range_alpha.data());
  this->a.sync();
  auto ga_dot = this->a.dot(this->sparse_array);
  double ref_dot = 0.;
  for (const auto &item : this->sparse_array) {
    ref_dot += this->range_alpha[item.first] * item.second;
  }
  {
    auto proxy = this->lock.scope();
    ASSERT_THAT(ga_dot, DoubleEq(ref_dot));
  }
  this->a.sync();
}

TYPED_TEST_P(DistrArrayCollectiveLinAlgF, times) {
  this->a.fill(this->alpha);
  this->b.fill(this->beta);
  TypeParam c{this->a};
  c.times(this->a, this->b);
  this->a.sync();
  auto from_ga_buffer_a = this->a.vec();
  auto from_ga_buffer_b = this->b.vec();
  auto from_ga_buffer_c = c.vec();
  auto ref_a = std::vector<double>(this->dim, this->alpha);
  auto ref_b = std::vector<double>(this->dim, this->beta);
  auto ref_c = std::vector<double>(this->dim, this->alpha * this->beta);
  {
    auto proxy = this->lock.scope();
    ASSERT_THAT(from_ga_buffer_a, Pointwise(DoubleEq(), ref_a));
    ASSERT_THAT(from_ga_buffer_b, Pointwise(DoubleEq(), ref_b));
    ASSERT_THAT(from_ga_buffer_c, Pointwise(DoubleEq(), ref_c));
  }
  this->a.sync();
}

// c[i] -= a[i]/(b[i]+shift)
TYPED_TEST_P(DistrArrayCollectiveLinAlgF, divide_append_negative) {
  double shift = 0.5;
  this->a.fill(this->alpha);
  this->b.fill(this->beta);
  TypeParam c{this->a};
  c.divide(this->a, this->b, shift, true, true);
  this->a.sync();
  auto from_ga_buffer_a = this->a.vec();
  auto from_ga_buffer_b = this->b.vec();
  auto from_ga_buffer_c = c.vec();
  auto ref_a = std::vector<double>(this->dim, this->alpha);
  auto ref_b = std::vector<double>(this->dim, this->beta);
  auto ref_c = std::vector<double>(this->dim, this->alpha - this->alpha / (this->beta + shift));
  {
    auto proxy = this->lock.scope();
    ASSERT_THAT(from_ga_buffer_a, Pointwise(DoubleEq(), ref_a));
    ASSERT_THAT(from_ga_buffer_b, Pointwise(DoubleEq(), ref_b));
    ASSERT_THAT(from_ga_buffer_c, Pointwise(DoubleEq(), ref_c));
  }
  this->a.sync();
}

// c[i] += a[i]/(b[i]+shift)
TYPED_TEST_P(DistrArrayCollectiveLinAlgF, divide_append_positive) {
  double shift = 0.5;
  this->a.fill(this->alpha);
  this->b.fill(this->beta);
  TypeParam c{this->a};
  c.divide(this->a, this->b, shift, true, false);
  this->a.sync();
  auto from_ga_buffer_a = this->a.vec();
  auto from_ga_buffer_b = this->b.vec();
  auto from_ga_buffer_c = c.vec();
  auto ref_a = std::vector<double>(this->dim, this->alpha);
  auto ref_b = std::vector<double>(this->dim, this->beta);
  auto ref_c = std::vector<double>(this->dim, this->alpha + this->alpha / (this->beta + shift));
  {
    auto proxy = this->lock.scope();
    ASSERT_THAT(from_ga_buffer_a, Pointwise(DoubleEq(), ref_a));
    ASSERT_THAT(from_ga_buffer_b, Pointwise(DoubleEq(), ref_b));
    ASSERT_THAT(from_ga_buffer_c, Pointwise(DoubleEq(), ref_c));
  }
  this->a.sync();
}

// c[i] = a[i]/(b[i]+shift)
TYPED_TEST_P(DistrArrayCollectiveLinAlgF, divide_overwrite_positive) {
  double shift = 0.5;
  this->a.fill(this->alpha);
  this->b.fill(this->beta);
  TypeParam c{this->a};
  c.divide(this->a, this->b, shift, false, false);
  this->a.sync();
  auto from_ga_buffer_a = this->a.vec();
  auto from_ga_buffer_b = this->b.vec();
  auto from_ga_buffer_c = c.vec();
  auto ref_a = std::vector<double>(this->dim, this->alpha);
  auto ref_b = std::vector<double>(this->dim, this->beta);
  auto ref_c = std::vector<double>(this->dim, this->alpha / (this->beta + shift));
  {
    auto proxy = this->lock.scope();
    ASSERT_THAT(from_ga_buffer_a, Pointwise(DoubleEq(), ref_a));
    ASSERT_THAT(from_ga_buffer_b, Pointwise(DoubleEq(), ref_b));
    ASSERT_THAT(from_ga_buffer_c, Pointwise(DoubleEq(), ref_c));
  }
  this->a.sync();
}

REGISTER_TYPED_TEST_SUITE_P(DistArrayBasicF, size, zero, fill);
REGISTER_TYPED_TEST_SUITE_P(DistArrayBasicRMAF, vec, get, put);
REGISTER_TYPED_TEST_SUITE_P(DistrArrayRangeRMAF, gather, scatter, scatter_acc, at);
REGISTER_TYPED_TEST_SUITE_P(DistrArrayRangeMinMaxF, min_loc_n, min_loc_n_reverse, max_n, min_abs_n, max_abs_n);
REGISTER_TYPED_TEST_SUITE_P(DistrArrayRangeLinAlgF, scal_double, add_double, sub_double, recip);
REGISTER_TYPED_TEST_SUITE_P(TestDistrArray, constructor, constructor_copy, constructor_copy_allocated,
                            copy_assignment_op_allocated, constructor_move, constructor_move_allocated,
                            move_assignment_op_allocated, select_max_dot);
REGISTER_TYPED_TEST_SUITE_P(DistrArrayCollectiveLinAlgF, add, sub, axpy, axpy_map, dot_array, dot_map, times,
                            divide_append_negative, divide_append_positive, divide_overwrite_positive);

#endif // LINEARALGEBRA_TEST_ARRAY_TESTDISTRARRAY_H
