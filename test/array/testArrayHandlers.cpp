#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "parallel_util.h"
#include <list>
#include <molpro/linalg/array/DistrArrayFile.h>
#include <molpro/linalg/array/DistrArrayMPI3.h>
#include <molpro/linalg/array/util.h>
#include <molpro/linalg/itsolv/ArrayHandlers.h>

using molpro::linalg::array::DistrArrayFile;
using molpro::linalg::array::DistrArrayMPI3;
using molpro::linalg::array::util::LockMPI3;
using molpro::linalg::itsolv::ArrayHandlers;
using molpro::linalg::test::mpi_comm;

using ::testing::ContainerEq;
using ::testing::DoubleEq;
using ::testing::Each;
using ::testing::Pointwise;

TEST(TestArrayHandlersIterable, constructor) {
  auto handlers = ArrayHandlers<std::vector<double>, std::vector<double>, std::map<size_t, double>>();
}

TEST(TestArrayHandlersDistrMPI3Sparse, constructor) {
  auto handlers = ArrayHandlers<DistrArrayMPI3, DistrArrayMPI3, std::map<size_t, double>>();
}

TEST(TestArrayHandlersDistrMPI3Sparse, axpy) {
  const size_t dim = 20;
  const double alpha = 2.;
  const double beta = 0.5;
  auto y = DistrArrayMPI3(dim, mpi_comm);
  y.fill(beta);
  auto y_ref = std::vector<double>(dim, beta);
  auto x = std::map<size_t, double>{{1, 1.0}, {3, 2.0}, {6, 3.0}, {11, 4.0}};
  for (const auto& el : x) {
    y_ref[el.first] += alpha * el.second;
  }
  auto handlers = ArrayHandlers<decltype(y), decltype(y), decltype(x)>();
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
  x.fill(alpha);
  x.sync();
  auto y = std::map<size_t, double>{{1, 1.0}, {3, 2.0}, {6, 3.0}, {11, 4.0}};
  auto handlers = ArrayHandlers<decltype(x), decltype(x), decltype(y)>();
  double ref = 0.5 + 0.5 * 2. + 0.5 * 3. + 0.5 * 4.;
  auto result = handlers.rp().dot(x, y);
  LockMPI3 lock{mpi_comm};
  auto l = lock.scope();
  ASSERT_DOUBLE_EQ(result, ref);
}

TEST(TestArrayHandlersIterableDDisk, default_handlers) {
  using molpro::linalg::array::ArrayHandlerDDiskIterable;
  using molpro::linalg::array::ArrayHandlerDefault;
  using molpro::linalg::array::ArrayHandlerIterableDDisk;
  using molpro::linalg::array::default_handler_t;
  using molpro::linalg::array::Span;
  using vec = std::vector<double>;
  static_assert(std::is_same_v<default_handler_t<vec, DistrArrayFile>, ArrayHandlerIterableDDisk<vec, DistrArrayFile>>);
  static_assert(std::is_same_v<default_handler_t<DistrArrayFile, vec>, ArrayHandlerDDiskIterable<DistrArrayFile, vec>>);
  static_assert(std::is_same_v<default_handler_t<Span<double>, DistrArrayFile>,
                               ArrayHandlerIterableDDisk<Span<double>, DistrArrayFile>>);
  static_assert(std::is_same_v<default_handler_t<DistrArrayFile, Span<double>>,
                               ArrayHandlerDDiskIterable<DistrArrayFile, Span<double>>>);
  // iterable arrays that are not contiguous are left to the default handler
  using list = std::list<double>;
  static_assert(std::is_same_v<default_handler_t<list, DistrArrayFile>, ArrayHandlerDefault<list, DistrArrayFile>>);
  static_assert(std::is_same_v<default_handler_t<DistrArrayFile, list>, ArrayHandlerDefault<DistrArrayFile, list>>);
}

TEST(TestArrayHandlersIterableDDisk, constructor) {
  auto handlers_rq = ArrayHandlers<std::vector<double>, DistrArrayFile, std::map<size_t, double>>();
  auto handlers_qr = ArrayHandlers<DistrArrayFile, std::vector<double>, std::map<size_t, double>>();
}

TEST(TestArrayHandlersIterableDDisk, rq_qr) {
  const size_t dim = 13;
  auto r = std::vector<double>(dim);
  std::iota(r.begin(), r.end(), 1.);
  auto handlers = ArrayHandlers<std::vector<double>, DistrArrayFile, std::map<size_t, double>>();
  // R is replicated, so the Q made from it lives on a single process
  auto q = handlers.qr().copy(r);
  EXPECT_EQ(q.communicator(), molpro::mpi::comm_self());
  EXPECT_DOUBLE_EQ(handlers.rq().dot(r, q), std::inner_product(r.begin(), r.end(), r.begin(), 0.));
  handlers.rq().axpy(-1., q, r);
  EXPECT_THAT(r, Each(DoubleEq(0.)));
}
