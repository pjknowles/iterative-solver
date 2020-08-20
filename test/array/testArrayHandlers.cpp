#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "parallel_util.h"
#include <molpro/linalg/array/DistrArrayMPI3.h>
#include <molpro/linalg/array/util.h>
#include <molpro/linalg/iterativesolver/ArrayHandlers.h>

using molpro::linalg::array::DistrArrayMPI3;
using molpro::linalg::array::util::LockMPI3;
using molpro::linalg::iterativesolver::ArrayHandlers;
using molpro::linalg::test::mpi_comm;

using ::testing::ContainerEq;
using ::testing::DoubleEq;
using ::testing::Each;
using ::testing::Pointwise;

TEST(TestArrayHandlersIterable, constructor) {
  auto handlers = ArrayHandlers<std::vector<double>, std::vector<double>, std::map<size_t, double>>::Builder{}.build();
}

TEST(TestArrayHandlersDistrMPI3Sparse, constructor) {
  auto handlers = ArrayHandlers<DistrArrayMPI3, DistrArrayMPI3, std::map<size_t, double>>::Builder{}.build();
}

TEST(TestArrayHandlersDistrMPI3Sparse, axpy) {
  const size_t dim = 20;
  const double alpha = 2.;
  const double beta = 0.5;
  auto y = DistrArrayMPI3(dim, mpi_comm);
  y.allocate_buffer();
  y.fill(beta);
  auto y_ref = std::vector<double>(dim, beta);
  auto x = std::map<size_t, double>{{1, 1.0}, {3, 2.0}, {6, 3.0}, {11, 4.0}};
  for (const auto& el : x) {
    y_ref[el.first] += alpha * el.second;
  }
  auto handlers = ArrayHandlers<decltype(y), decltype(y), decltype(x)>::Builder{}.build();
  handlers.rp().axpy(alpha, x, y);
  y.sync();
  LockMPI3 lock{mpi_comm};
  auto l = lock.scope();
  auto y_result = y.vec();
  ASSERT_THAT(y_result, Pointwise(DoubleEq(), y_ref));
}

TEST(TestArrayHandlersDistrMPI3Sparse, dot) {
  const size_t dim = 20;
  const double alpha = 0.5;
  auto x = DistrArrayMPI3(dim, mpi_comm);
  x.allocate_buffer();
  x.fill(alpha);
  x.sync();
  auto y = std::map<size_t, double>{{1, 1.0}, {3, 2.0}, {6, 3.0}, {11, 4.0}};
  auto handlers = ArrayHandlers<decltype(x), decltype(x), decltype(y)>::Builder{}.build();
  double ref = 0.5 + 0.5 * 2. + 0.5 * 3. + 0.5 * 4.;
  auto result = handlers.rp().dot(x, y);
  LockMPI3 lock{mpi_comm};
  auto l = lock.scope();
  ASSERT_DOUBLE_EQ(result, ref);
}
