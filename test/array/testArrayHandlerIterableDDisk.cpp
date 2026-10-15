#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <numeric>
#include <type_traits>
#include <vector>

#include "parallel_util.h"
#include <molpro/linalg/array/ArrayHandlerDDiskIterable.h>
#include <molpro/linalg/array/ArrayHandlerIterable.h>
#include <molpro/linalg/array/ArrayHandlerIterableDDisk.h>
#include <molpro/linalg/array/DistrArrayFile.h>
#include <molpro/linalg/array/DistrArrayMPI3.h>
#include <molpro/linalg/array/Span.h>
#include <molpro/linalg/array/util/gemm.h>
#include <molpro/linalg/itsolv/wrap.h>
#include <molpro/linalg/options.h>
#include <molpro/mpi.h>

using molpro::linalg::array::ArrayHandlerDDiskIterable;
using molpro::linalg::array::ArrayHandlerIterable;
using molpro::linalg::array::ArrayHandlerIterableDDisk;
using molpro::linalg::array::DistrArrayFile;
using molpro::linalg::array::DistrArrayMPI3;
using molpro::linalg::array::Span;
using molpro::linalg::itsolv::cwrap;
using molpro::linalg::itsolv::wrap;
using molpro::linalg::itsolv::subspace::Matrix;
using molpro::linalg::test::mpi_comm;

namespace {
//! The full content of a distributed array, assembled from the local sections of all processes
template <class Q>
std::vector<double> full(const Q& q) {
  auto result = std::vector<double>(q.size(), 0.);
  const auto [lo, hi] = molpro::linalg::array::util::local_range(q);
  if (hi > lo)
    q.get(lo, hi, result.data() + lo);
  MPI_Allreduce(MPI_IN_PLACE, result.data(), result.size(), MPI_DOUBLE, MPI_SUM, q.communicator());
  return result;
}

//! Sets the content of a distributed array from the full array v, each process putting its local section
template <class Q>
void set(Q& q, const std::vector<double>& v) {
  const auto [lo, hi] = molpro::linalg::array::util::local_range(q);
  if (hi > lo)
    q.put(lo, hi, v.data() + lo);
}

void expect_near(const std::vector<double>& result, const std::vector<double>& reference) {
  ASSERT_EQ(result.size(), reference.size());
  for (size_t i = 0; i < result.size(); ++i)
    EXPECT_NEAR(result[i], reference[i], 1e-11 * (1 + std::abs(reference[i]))) << "element " << i;
}

void expect_near(const Matrix<double>& result, const Matrix<double>& reference) {
  ASSERT_EQ(result.rows(), reference.rows());
  ASSERT_EQ(result.cols(), reference.cols());
  for (size_t i = 0; i < result.rows(); ++i)
    for (size_t j = 0; j < result.cols(); ++j)
      EXPECT_NEAR(result(i, j), reference(i, j), 1e-11 * (1 + std::abs(reference(i, j)))) << i << ", " << j;
}
} // namespace

//! A test configuration: the contiguous array type, and whether the distributed arrays are on comm_self
template <class R, bool on_comm_self>
struct Config {
  using Rtype = R;
  static constexpr bool self = on_comm_self;
};

/*!
 * Contiguous arrays rr and distributed arrays qq (DistrArrayFile), with plain vectors holding the same values as
 * references. std::vector arrays are separate, so not evenly spaced, while Span arrays are views into one buffer, and
 * so are evenly spaced, which exercises both the dgemv and the dgemm paths of gemm.
 */
template <class C>
class ArrayHandlerIterableDDiskF : public ::testing::Test {
public:
  using R = typename C::Rtype;
  using Q = DistrArrayFile;
  static constexpr size_t dim = 47; // not a multiple of the page sizes
  static constexpr size_t nr = 3;
  static constexpr size_t nq = 4;

  ArrayHandlerIterableDDiskF() : comm(C::self ? molpro::mpi::comm_self() : mpi_comm) {
    molpro::linalg::set_options(molpro::Options("ITERATIVE-SOLVER", "GEMM_PAGESIZE=5"));
    for (size_t i = 0; i < nr; ++i) {
      ref_r.emplace_back(dim);
      for (size_t k = 0; k < dim; ++k)
        ref_r.back()[k] = std::sin(1. + i + 0.37 * k);
    }
    for (size_t j = 0; j < nq; ++j) {
      ref_q.emplace_back(dim);
      for (size_t k = 0; k < dim; ++k)
        ref_q.back()[k] = std::cos(0.5 + 2. * j + 0.23 * k) - 0.1 * j;
    }
    if constexpr (std::is_same_v<R, std::vector<double>>) {
      rr = ref_r;
    } else {
      buffer.resize(nr * dim);
      for (size_t i = 0; i < nr; ++i) {
        std::copy(ref_r[i].begin(), ref_r[i].end(), buffer.begin() + i * dim);
        rr.emplace_back(buffer.data() + i * dim, dim);
      }
    }
    for (size_t j = 0; j < nq; ++j) {
      qq.emplace_back(dim, comm);
      qq.back().set_buffer_size(7);
      set(qq.back(), ref_q[j]);
    }
  }
  ~ArrayHandlerIterableDDiskF() override { molpro::linalg::set_options(molpro::Options("ITERATIVE-SOLVER", "")); }

  static std::vector<double> to_vector(const R& r) { return {std::data(r), std::data(r) + std::size(r)}; }

  MPI_Comm comm;
  std::vector<std::vector<double>> ref_r, ref_q;
  std::vector<double> buffer;
  std::vector<R> rr;
  std::vector<Q> qq;
  ArrayHandlerIterableDDisk<R, Q> rq;
  ArrayHandlerDDiskIterable<Q, R> qr;
  ArrayHandlerIterable<std::vector<double>> ref;
};

using Configs = ::testing::Types<Config<std::vector<double>, true>, Config<std::vector<double>, false>,
                                 Config<Span<double>, true>, Config<Span<double>, false>>;
TYPED_TEST_SUITE(ArrayHandlerIterableDDiskF, Configs);

TYPED_TEST(ArrayHandlerIterableDDiskF, copy) {
  this->rq.copy(this->rr[0], this->qq[1]);
  expect_near(this->to_vector(this->rr[0]), this->ref_q[1]);
  this->qr.copy(this->qq[2], this->rr[1]);
  expect_near(full(this->qq[2]), this->ref_r[1]);
}

TYPED_TEST(ArrayHandlerIterableDDiskF, round_trip) {
  auto q = this->qr.copy(this->rr[2]);
  this->rq.fill(0, this->rr[0]);
  this->rq.copy(this->rr[0], q);
  expect_near(this->to_vector(this->rr[0]), this->ref_r[2]);
}

TYPED_TEST(ArrayHandlerIterableDDiskF, construct_copy) {
  if constexpr (std::is_same_v<typename TestFixture::R, std::vector<double>>) {
    auto r = this->rq.copy(this->qq[3]);
    expect_near(r, this->ref_q[3]);
  }
  // R is not distributed, so neither is a Q made from it
  auto q = this->qr.copy(this->rr[1]);
  int comparison;
  MPI_Comm_compare(q.communicator(), molpro::mpi::comm_self(), &comparison);
  EXPECT_TRUE(comparison == MPI_IDENT or comparison == MPI_CONGRUENT);
  EXPECT_EQ(q.size(), this->dim);
  expect_near(full(q), this->ref_r[1]);
}

TYPED_TEST(ArrayHandlerIterableDDiskF, construct_copy_directory) {
  int rank;
  MPI_Comm_rank(mpi_comm, &rank);
  const auto directory = std::filesystem::absolute("ArrayHandlerIterableDDisk-directory-" + std::to_string(rank));
  std::filesystem::create_directories(directory);
  {
    auto handler = ArrayHandlerDDiskIterable<DistrArrayFile, typename TestFixture::R>(directory.string());
    auto q = handler.copy(this->rr[0]);
    expect_near(full(q), this->ref_r[0]);
  }
  std::filesystem::remove_all(directory);
}

TYPED_TEST(ArrayHandlerIterableDDiskF, construct_copy_function) {
  int calls = 0;
  auto handler =
      ArrayHandlerDDiskIterable<DistrArrayFile, typename TestFixture::R>([&calls](const typename TestFixture::R& r) {
        ++calls;
        auto q = DistrArrayFile(std::size(r), molpro::mpi::comm_self());
        q.fill(42);
        return q;
      });
  auto q = handler.copy(this->rr[0]);
  EXPECT_EQ(calls, 1);
  EXPECT_THAT(full(q), ::testing::Each(::testing::DoubleEq(42)));
  EXPECT_EQ(handler.counter().copy, 1);
}

TYPED_TEST(ArrayHandlerIterableDDiskF, scal_fill) {
  this->rq.scal(2.5, this->rr[0]);
  this->ref.scal(2.5, this->ref_r[0]);
  expect_near(this->to_vector(this->rr[0]), this->ref_r[0]);
  this->rq.fill(-1.5, this->rr[1]);
  expect_near(this->to_vector(this->rr[1]), std::vector<double>(this->dim, -1.5));
  this->qr.scal(-0.5, this->qq[0]);
  this->ref.scal(-0.5, this->ref_q[0]);
  expect_near(full(this->qq[0]), this->ref_q[0]);
  this->qr.fill(3., this->qq[1]);
  expect_near(full(this->qq[1]), std::vector<double>(this->dim, 3.));
}

TYPED_TEST(ArrayHandlerIterableDDiskF, axpy) {
  this->rq.axpy(-1.25, this->qq[2], this->rr[1]);
  this->ref.axpy(-1.25, this->ref_q[2], this->ref_r[1]);
  expect_near(this->to_vector(this->rr[1]), this->ref_r[1]);
  this->qr.axpy(0.75, this->rr[2], this->qq[3]);
  this->ref.axpy(0.75, this->ref_r[2], this->ref_q[3]);
  expect_near(full(this->qq[3]), this->ref_q[3]);
}

TYPED_TEST(ArrayHandlerIterableDDiskF, dot) {
  const auto reference = this->ref.dot(this->ref_r[2], this->ref_q[1]);
  EXPECT_NEAR(this->rq.dot(this->rr[2], this->qq[1]), reference, 1e-12);
  EXPECT_NEAR(this->qr.dot(this->qq[1], this->rr[2]), reference, 1e-12);
}

TYPED_TEST(ArrayHandlerIterableDDiskF, gemm_inner) {
  const auto reference = this->ref.gemm_inner(cwrap(this->ref_r), cwrap(this->ref_q));
  const auto result_rq = this->rq.gemm_inner(cwrap(this->rr), cwrap(this->qq));
  expect_near(result_rq, reference);
  const auto result_qr = this->qr.gemm_inner(cwrap(this->qq), cwrap(this->rr));
  expect_near(result_qr, molpro::linalg::array::util::conjugate_transpose(reference));
  EXPECT_EQ(this->rq.counter().gemm_inner, 1);
  EXPECT_EQ(this->qr.counter().gemm_inner, 1);
}

TYPED_TEST(ArrayHandlerIterableDDiskF, gemm_outer) {
  // rr += alpha^T qq
  auto alpha_qr = Matrix<double>({this->nq, this->nr});
  for (size_t i = 0; i < alpha_qr.rows(); ++i)
    for (size_t j = 0; j < alpha_qr.cols(); ++j)
      alpha_qr(i, j) = 0.5 + i - 0.25 * j;
  this->rq.gemm_outer(alpha_qr, cwrap(this->qq), wrap(this->rr));
  this->ref.gemm_outer(alpha_qr, cwrap(this->ref_q), wrap(this->ref_r));
  for (size_t i = 0; i < this->nr; ++i)
    expect_near(this->to_vector(this->rr[i]), this->ref_r[i]);
  // qq += alpha^T rr
  auto alpha_rq = Matrix<double>({this->nr, this->nq});
  for (size_t i = 0; i < alpha_rq.rows(); ++i)
    for (size_t j = 0; j < alpha_rq.cols(); ++j)
      alpha_rq(i, j) = -1. + 0.3 * i + j;
  this->qr.gemm_outer(alpha_rq, cwrap(this->rr), wrap(this->qq));
  this->ref.gemm_outer(alpha_rq, cwrap(this->ref_r), wrap(this->ref_q));
  for (size_t j = 0; j < this->nq; ++j)
    expect_near(full(this->qq[j]), this->ref_q[j]);
}

TYPED_TEST(ArrayHandlerIterableDDiskF, gemm_empty) {
  using R = typename TestFixture::R;
  using Q = typename TestFixture::Q;
  auto no_r = std::vector<R>();
  auto no_q = std::vector<Q>();
  EXPECT_EQ(this->rq.gemm_inner(cwrap(no_r), cwrap(this->qq)).size(), 0);
  EXPECT_EQ(this->qr.gemm_inner(cwrap(no_q), cwrap(this->rr)).size(), 0);
  this->rq.gemm_outer(Matrix<double>({0, this->nr}), cwrap(no_q), wrap(this->rr));
  this->qr.gemm_outer(Matrix<double>({0, this->nq}), cwrap(no_r), wrap(this->qq));
  for (size_t i = 0; i < this->nr; ++i)
    expect_near(this->to_vector(this->rr[i]), this->ref_r[i]);
  for (size_t j = 0; j < this->nq; ++j)
    expect_near(full(this->qq[j]), this->ref_q[j]);
}

// An empty contiguous array need not have any memory, i.e. its data() may be nullptr
TYPED_TEST(ArrayHandlerIterableDDiskF, empty_arrays) {
  using R = typename TestFixture::R;
  using Q = typename TestFixture::Q;
  auto rr = std::vector<R>(2);
  auto qq = std::vector<Q>();
  qq.emplace_back(0, this->comm);
  qq.emplace_back(0, this->comm);
  this->rq.copy(rr[0], qq[0]);
  this->qr.copy(qq[0], rr[0]);
  this->rq.axpy(2., qq[0], rr[0]);
  this->qr.axpy(2., rr[0], qq[0]);
  EXPECT_EQ(this->rq.dot(rr[0], qq[0]), 0.);
  EXPECT_EQ(this->qr.dot(qq[0], rr[0]), 0.);
  EXPECT_TRUE(this->rq.select_max_dot(0, rr[0], qq[0]).empty());
  EXPECT_TRUE(this->qr.select_max_dot(0, qq[0], rr[0]).empty());
  EXPECT_EQ(this->qr.copy(rr[0]).size(), 0);
  if constexpr (std::is_same_v<R, std::vector<double>>)
    EXPECT_TRUE(this->rq.copy(qq[0]).empty());
  const auto inner = this->rq.gemm_inner(cwrap(rr), cwrap(qq));
  ASSERT_EQ(inner.rows(), 2);
  ASSERT_EQ(inner.cols(), 2);
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      EXPECT_EQ(inner(i, j), 0.);
  EXPECT_EQ(this->qr.gemm_inner(cwrap(qq), cwrap(rr)).size(), 4);
  auto alphas = Matrix<double>({2, 2});
  alphas.fill(1.);
  this->rq.gemm_outer(alphas, cwrap(qq), wrap(rr));
  this->qr.gemm_outer(alphas, cwrap(rr), wrap(qq));
}

TYPED_TEST(ArrayHandlerIterableDDiskF, select_max_dot) {
  const size_t n = 5;
  const auto reference = this->ref.select_max_dot(n, this->ref_r[0], this->ref_q[0]);
  const auto result_rq = this->rq.select_max_dot(n, this->rr[0], this->qq[0]);
  const auto result_qr = this->qr.select_max_dot(n, this->qq[0], this->rr[0]);
  ASSERT_EQ(result_rq.size(), n);
  for (const auto& [i, v] : reference) {
    ASSERT_EQ(result_rq.count(i), 1) << "index " << i;
    EXPECT_NEAR(result_rq.at(i), v, 1e-12);
  }
  EXPECT_EQ(result_rq, result_qr);
}

TYPED_TEST(ArrayHandlerIterableDDiskF, select) {
  const size_t n = 4;
  for (bool max : {false, true})
    for (bool ignore_sign : {false, true}) {
      const auto reference_r = this->ref.select(n, this->ref_r[1], max, ignore_sign);
      EXPECT_EQ(this->rq.select(n, this->rr[1], max, ignore_sign), reference_r);
      const auto reference_q = this->ref.select(n, this->ref_q[1], max, ignore_sign);
      const auto result_q = this->qr.select(n, this->qq[1], max, ignore_sign);
      ASSERT_EQ(result_q.size(), n);
      for (const auto& [i, v] : reference_q) {
        ASSERT_EQ(result_q.count(i), 1) << "index " << i << " max " << max << " ignore_sign " << ignore_sign;
        EXPECT_NEAR(result_q.at(i), v, 1e-12);
      }
    }
}

TYPED_TEST(ArrayHandlerIterableDDiskF, lazy_dot_and_axpy) {
  auto result = std::vector<double>(this->nr * this->nq, 0.);
  {
    auto h = this->rq.lazy_handle();
    for (size_t i = 0; i < this->nr; ++i)
      for (size_t j = 0; j < this->nq; ++j)
        h.dot(this->rr[i], this->qq[j], result[i * this->nq + j]);
  }
  for (size_t i = 0; i < this->nr; ++i)
    for (size_t j = 0; j < this->nq; ++j)
      EXPECT_NEAR(result[i * this->nq + j], this->ref.dot(this->ref_r[i], this->ref_q[j]), 1e-12);
  {
    auto h = this->qr.lazy_handle();
    h.axpy(2., this->rr[0], this->qq[1]);
    h.axpy(-1., this->rr[2], this->qq[1]);
  }
  this->ref.axpy(2., this->ref_r[0], this->ref_q[1]);
  this->ref.axpy(-1., this->ref_r[2], this->ref_q[1]);
  expect_near(full(this->qq[1]), this->ref_q[1]);
}

// The handlers work with any distributed array; here with one held in memory
TEST(ArrayHandlerIterableDDisk, distributed_in_memory) {
  const size_t dim = 23;
  auto r = std::vector<double>(dim), x = std::vector<double>(dim);
  std::iota(r.begin(), r.end(), 1.);
  for (size_t k = 0; k < dim; ++k)
    x[k] = 0.5 - 0.1 * k;
  auto q = DistrArrayMPI3(dim, mpi_comm);
  set(q, x);
  auto rq = ArrayHandlerIterableDDisk<std::vector<double>, DistrArrayMPI3>();
  auto qr = ArrayHandlerDDiskIterable<DistrArrayMPI3, std::vector<double>>([](const std::vector<double>& source) {
    auto result = DistrArrayMPI3(source.size(), mpi_comm);
    result.copy(Span<const double>(source.data(), source.size()));
    return result;
  });
  auto ref = ArrayHandlerIterable<std::vector<double>>();
  EXPECT_NEAR(rq.dot(r, q), ref.dot(r, x), 1e-12);
  rq.axpy(2., q, r);
  auto ref_r = std::vector<double>(dim);
  std::iota(ref_r.begin(), ref_r.end(), 1.);
  ref.axpy(2., x, ref_r);
  expect_near(r, ref_r);
  auto qq = std::vector<DistrArrayMPI3>{qr.copy(r)};
  expect_near(full(qq[0]), r);
  auto alpha = Matrix<double>({1, 1});
  alpha(0, 0) = -3.;
  auto xx = std::vector<std::vector<double>>{x};
  qr.gemm_outer(alpha, cwrap(xx), wrap(qq));
  ref.axpy(-3., x, ref_r);
  expect_near(full(qq[0]), ref_r);
  qq[0].sync();
}
