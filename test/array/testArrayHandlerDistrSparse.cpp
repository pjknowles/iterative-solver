#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "parallel_util.h"
#include <molpro/linalg/array/ArrayHandlerDistrSparse.h>
#include <molpro/linalg/array/DistrArrayMPI3.h>
#include <molpro/linalg/array/DistrArraySpan.h>
#include <molpro/linalg/array/DistrArrayFile.h>
#include <molpro/linalg/array/util.h>
#include <molpro/linalg/itsolv/ArrayHandlers.h>

using molpro::linalg::array::ArrayHandlerDistrSparse;
using molpro::linalg::array::DistrArrayMPI3;
using molpro::linalg::itsolv::ArrayHandlers;
using molpro::linalg::array::DistrArraySpan;
using molpro::linalg::array::DistrArrayFile;
using molpro::linalg::array::util::LockMPI3;
using molpro::linalg::test::mpi_comm;

using ::testing::ContainerEq;
using ::testing::DoubleEq;
using ::testing::Each;
using ::testing::Pointwise;

TEST(TestArrayHandlerDistrSparse, constructor) {
  auto handler = ArrayHandlerDistrSparse<DistrArrayMPI3, std::map<size_t, double>>{};
}

//TEST(TestArrayHandlerDistrSpanSparse, constructor) {
//  //auto handler = ArrayHandlerDistrSparse<DistrArraySpan, std::map<size_t, double>>{};
//  const std::shared_ptr<ArrayHandlers<DistrArraySpan, DistrArrayFile, std::map<size_t, double>>>& handlers =
//                                                   std::make_shared<ArrayHandlers<DistrArraySpan,
//                                                                    DistrArrayFile, std::map<size_t, double>>>();
//}

TEST(TestArrayHandlerDistrSparse, axpy) {
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
  auto handler = ArrayHandlerDistrSparse<decltype(y), decltype(x)>{};
  handler.axpy(alpha, x, y);
  y.sync();
  LockMPI3 lock{mpi_comm};
  auto l = lock.scope();
  auto y_result = y.vec();
  ASSERT_THAT(y_result, Pointwise(DoubleEq(), y_ref));
}

TEST(TestArrayHandlerDistrSparse, dot) {
  const size_t dim = 20;
  const double alpha = 0.5;
  auto x = DistrArrayMPI3(dim, mpi_comm);
  x.allocate_buffer();
  x.fill(alpha);
  x.sync();
  auto y = std::map<size_t, double>{{1, 1.0}, {3, 2.0}, {6, 3.0}, {11, 4.0}};
  auto handler = ArrayHandlerDistrSparse<decltype(x), decltype(y)>{};
  double ref = 0.5 + 0.5 * 2. + 0.5 * 3. + 0.5 * 4.;
  auto result = handler.dot(x, y);
  LockMPI3 lock{mpi_comm};
  auto l = lock.scope();
  ASSERT_DOUBLE_EQ(result, ref);
}
