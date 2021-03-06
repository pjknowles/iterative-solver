#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <numeric>
#include <utility>

#include "parallel_util.h"
#include <molpro/linalg/array/ArrayHandlerDistr.h>
#include <molpro/linalg/array/DistrArrayMPI3.h>
#include <molpro/linalg/array/util.h>

using molpro::linalg::array::ArrayHandlerDistr;
using molpro::linalg::array::DistrArrayMPI3;
using molpro::linalg::array::util::LockMPI3;
using molpro::linalg::test::mpi_comm;

using ::testing::ContainerEq;
using ::testing::Each;

TEST(TestArrayHandlerDistr, lazy_dot) {
  LockMPI3 lock{mpi_comm};
  ArrayHandlerDistr<DistrArrayMPI3> handler{};
  using value_type = double;
  static const int N = 3;
  static const int dim = 11;
  std::vector<DistrArrayMPI3> xx{N, DistrArrayMPI3(dim, mpi_comm)};
  std::vector<DistrArrayMPI3> yy{N, DistrArrayMPI3(dim, mpi_comm)};
  auto result = std::vector<value_type>(N * N);
  auto ref_dot = std::vector<value_type>(N * N);
  auto vec = std::vector<double>(dim);
  for (size_t i = 0; i < N; ++i) {
    auto x = xx[i].local_buffer();
    std::iota(begin(*x), end(*x), (value_type)i + x->start());
    auto y = yy[i].local_buffer();
    std::iota(begin(*y), end(*y), (value_type)100 * i + y->start());
  }
  {
    auto h = handler.lazy_handle();
    for (size_t i = 0, ij = 0; i < N; ++i) {
      for (size_t j = 0; j < N; ++j, ++ij) {
        ref_dot[ij] = handler.dot(xx[i], yy[j]);
        h.dot(xx[i], yy[j], result[ij]);
      }
    }
    {
      auto l = lock.scope();
      EXPECT_FALSE(h.invalid());
      EXPECT_THAT(result, Each(0));
    }
  }
  {
    auto l = lock.scope();
    EXPECT_THAT(result, ContainerEq(ref_dot));
  }
}
