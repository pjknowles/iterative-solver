#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <vector>

#include "parallel_util.h"

#include <molpro/linalg/array/ArrayHandlerDDisk.h>
#include <molpro/linalg/array/ArrayHandlerDDiskDistr.h>
#include <molpro/linalg/array/ArrayHandlerDDiskSparse.h>
#include <molpro/linalg/array/ArrayHandlerDistr.h>
#include <molpro/linalg/array/ArrayHandlerDistrDDisk.h>
#include <molpro/linalg/array/ArrayHandlerDistrSparse.h>
#include <molpro/linalg/array/ArrayHandlerIterable.h>
#include <molpro/linalg/array/ArrayHandlerIterableSparse.h>
#include <molpro/linalg/array/ArrayHandlerSparse.h>
#include <molpro/linalg/array/DistrArrayDisk.h>
#include <molpro/linalg/array/DistrArrayFile.h>
#include <molpro/linalg/array/DistrArraySpan.h>
#ifdef LINEARALGEBRA_ARRAY_MPI3
#include <molpro/linalg/array/DistrArrayMPI3.h>
#endif
#include <molpro/linalg/array/Span.h>
#include <molpro/linalg/array/util.h>
#include <molpro/linalg/array/util/Distribution.h>
#include <molpro/linalg/itsolv/subspace/Matrix.h>
#include <molpro/linalg/itsolv/wrap.h>

using molpro::linalg::array::ArrayHandlerDDisk;
using molpro::linalg::array::ArrayHandlerDDiskDistr;
using molpro::linalg::array::ArrayHandlerDDiskSparse;
using molpro::linalg::array::ArrayHandlerDistr;
using molpro::linalg::array::ArrayHandlerDistrDDisk;
using molpro::linalg::array::ArrayHandlerDistrSparse;
using molpro::linalg::array::ArrayHandlerIterable;
using molpro::linalg::array::ArrayHandlerIterableSparse;
using molpro::linalg::array::ArrayHandlerSparse;
using molpro::linalg::array::DistrArrayDisk;
using molpro::linalg::array::DistrArrayFile;
using molpro::linalg::array::DistrArraySpan;
#ifdef LINEARALGEBRA_ARRAY_MPI3
using molpro::linalg::array::DistrArrayMPI3;
#endif
using molpro::linalg::array::Span;
using molpro::linalg::array::util::LockMPI3;
using molpro::linalg::array::util::make_distribution_spread_remainder;
using molpro::linalg::array::util::ScopeLock;
using molpro::linalg::itsolv::cwrap;
using molpro::linalg::itsolv::wrap;
using molpro::linalg::itsolv::subspace::Matrix;

using ::testing::ContainerEq;
using ::testing::DoubleEq;
using ::testing::Each;
using ::testing::Pointwise;

TEST(TestGemm, distr_inner) {
  using molpro::mpi::comm_global;

  auto handler = ArrayHandlerDistr<DistrArraySpan, DistrArraySpan>{};
  size_t n = 10;
  size_t dim = 10;
  std::vector<std::vector<double>> vx(n, std::vector<double>(dim)), vy(n, std::vector<double>(dim));
  std::vector<DistrArraySpan> cx, cy;
  cx.reserve(n);
  cy.reserve(n);
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm_global(), &mpi_rank);
  MPI_Comm_size(comm_global(), &mpi_size);
  for (size_t i = 0; i < n; i++) {
    std::iota(vx[i].begin(), vx[i].end(), i + 0.5);
    std::iota(vy[i].begin(), vy[i].end(), i + 0.5);
    auto crange = make_distribution_spread_remainder<size_t>(dim, mpi_size).range(mpi_rank);
    auto clength = crange.second - crange.first;
    cx.emplace_back(dim, Span<DistrArraySpan::value_type>(&vx[i][crange.first], clength), comm_global());
    cy.emplace_back(dim, Span<DistrArraySpan::value_type>(&vy[i][crange.first], clength), comm_global());
  }
  std::vector<double> vref(n * n), vgemm(n * n);
  std::pair<size_t, size_t> mat_dim = std::make_pair(n, n);
  Matrix<double> gemm_dot(vref, mat_dim);
  Matrix<double> ref_dot(vgemm, mat_dim);
  gemm_dot = handler.gemm_inner(cwrap(cx), cwrap(cy));
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      ref_dot(i, j) = handler.dot(cx[i], cy[j]);
    }
  }
  EXPECT_THAT(vgemm, Pointwise(DoubleEq(), vref));
}

TEST(TestGemm, distrarrayfile_inner) {
  using molpro::mpi::comm_global;

  auto handler = ArrayHandlerDistr<DistrArraySpan, DistrArrayFile>{};
  size_t n = 250;
  size_t dim = 250;
  std::vector<std::vector<double>> vx(n, std::vector<double>(dim)), vy(n, std::vector<double>(dim));
  std::vector<DistrArraySpan> cx;
  std::vector<DistrArrayFile> cy;
  cx.reserve(n);
  cy.reserve(n);

  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm_global(), &mpi_rank);
  MPI_Comm_size(comm_global(), &mpi_size);
  for (size_t i = 0; i < n; i++) {
    std::iota(vx[i].begin(), vx[i].end(), i + 0.5);
    std::iota(vy[i].begin(), vy[i].end(), i + 0.5);
    auto crange = make_distribution_spread_remainder<size_t>(dim, mpi_size).range(mpi_rank);
    auto clength = crange.second - crange.first;
    cx.emplace_back(dim, Span<DistrArraySpan::value_type>(&vx[i][crange.first], clength), comm_global());
    cy.emplace_back(dim, comm_global());
    cy.back().put(crange.first, crange.second, vy[i].data());
  }
  std::vector<double> vref, vgemm;
  Matrix<double> gemm_dot({n, n});
  gemm_dot = handler.gemm_inner(cwrap(cx), cwrap(cy));
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      vgemm.push_back(gemm_dot(i, j));
      vref.push_back(handler.dot(cx[i], cy[j]));
    }
  }
  EXPECT_THAT(vgemm, Pointwise(DoubleEq(), vref));
}

TEST(TestGemm, ddiskdistr_inner) {
  using molpro::mpi::comm_global;

  auto handler = ArrayHandlerDDiskDistr<DistrArrayFile, DistrArraySpan>{};
  size_t n = 10;
  size_t dim = 10;
  std::vector<std::vector<double>> vx(n, std::vector<double>(dim)), vy(n, std::vector<double>(dim));
  std::vector<DistrArrayFile> cx;
  std::vector<DistrArraySpan> cy;
  cx.reserve(n);
  cy.reserve(n);
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm_global(), &mpi_rank);
  MPI_Comm_size(comm_global(), &mpi_size);
  for (size_t i = 0; i < n; i++) {
    std::iota(vx[i].begin(), vx[i].end(), i + 0.5);
    std::iota(vy[i].begin(), vy[i].end(), i + 0.5);
    cx.emplace_back(dim);
    auto crange = cx.back().distribution().range(mpi_rank);
    auto clength = crange.second - crange.first;
    cx.back().put(crange.first, crange.second, &(*(vx[i].cbegin() + crange.first)));
    cy.emplace_back(dim, Span<DistrArraySpan::value_type>(&vy[i][crange.first], clength), comm_global());
  }

  std::vector<double> vref(n * n), vgemm(n * n);
  std::pair<size_t, size_t> mat_dim = std::make_pair(n, n);
  Matrix<double> gemm_dot(vref, mat_dim);
  Matrix<double> ref_dot(vgemm, mat_dim);
//  EXPECT_THROW(gemm_dot = handler.gemm_inner(cwrap(cx), cwrap(cy)), std::runtime_error);
  gemm_dot = handler.gemm_inner(cwrap(cx), cwrap(cy));
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      ref_dot(i, j) = handler.dot(cx[i], cy[j]);
    }
  }
  EXPECT_THAT(vgemm, Pointwise(DoubleEq(), vref));
}

TEST(TestGemm, distrddisk_inner) {
  using molpro::mpi::comm_global;

  auto handler = ArrayHandlerDistrDDisk<DistrArraySpan, DistrArrayFile>{};
  size_t n = 10;
  size_t dim = 10;
  std::vector<std::vector<double>> vx(n, std::vector<double>(dim)), vy(n, std::vector<double>(dim));
  std::vector<DistrArrayFile> cx;
  std::vector<DistrArraySpan> cy;
  cx.reserve(n);
  cy.reserve(n);
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm_global(), &mpi_rank);
  MPI_Comm_size(comm_global(), &mpi_size);
  for (size_t i = 0; i < n; i++) {
    std::iota(vx[i].begin(), vx[i].end(), i + 0.5);
    std::iota(vy[i].begin(), vy[i].end(), i + 0.5);
    cx.emplace_back(dim);
    auto crange = cx.back().distribution().range(mpi_rank);
    auto clength = crange.second - crange.first;
    cx.back().put(crange.first, crange.second, &(*(vx[i].cbegin() + crange.first)));
    cy.emplace_back(dim, Span<DistrArraySpan::value_type>(&vy[i][crange.first], clength), comm_global());
  }
  std::vector<double> vref(n * n), vgemm(n * n);
  std::pair<size_t, size_t> mat_dim = std::make_pair(n, n);
  Matrix<double> gemm_dot(vref, mat_dim);
  Matrix<double> ref_dot(vgemm, mat_dim);
  gemm_dot = handler.gemm_inner(cwrap(cy), cwrap(cx));
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      ref_dot(i, j) = handler.dot(cy[i], cx[j]);
    }
  }
  EXPECT_THAT(vgemm, Pointwise(DoubleEq(), vref));
}

TEST(TestGemm, ddisk_inner) {
  using molpro::mpi::comm_global;

  auto handler = ArrayHandlerDDisk<DistrArrayFile, DistrArrayFile>{};
  size_t n = 10;
  size_t dim = 10;
  std::vector<std::vector<double>> vx(n, std::vector<double>(dim)), vy(n, std::vector<double>(dim));
  std::vector<DistrArrayFile> cx;
  std::vector<DistrArrayFile> cy;
  cx.reserve(n);
  cy.reserve(n);
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm_global(), &mpi_rank);
  MPI_Comm_size(comm_global(), &mpi_size);
  for (size_t i = 0; i < n; i++) {
    std::iota(vx[i].begin(), vx[i].end(), i + 0.5);
    std::iota(vy[i].begin(), vy[i].end(), i + 0.5);
    cx.emplace_back(dim);
    cy.emplace_back(dim);
    auto crange = make_distribution_spread_remainder<size_t>(dim, mpi_size).range(mpi_rank);
    cx.back().put(crange.first, crange.second, &(*(vx[i].cbegin() + crange.first)));
    cy.back().put(crange.first, crange.second, &(*(vy[i].cbegin() + crange.first)));
  }
  std::vector<double> vref(n * n), vgemm(n * n);
  std::pair<size_t, size_t> mat_dim = std::make_pair(n, n);
  Matrix<double> gemm_dot(vref, mat_dim);
  Matrix<double> ref_dot(vgemm, mat_dim);
  gemm_dot = handler.gemm_inner(cwrap(cx), cwrap(cy));
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      ref_dot(i, j) = handler.dot(cx[i], cy[j]);
    }
  }
  EXPECT_THAT(vgemm, Pointwise(DoubleEq(), vref));
}

TEST(TestGemm, distrsparse_inner) {
  using molpro::mpi::comm_global;

  auto handler = ArrayHandlerDistrSparse<DistrArraySpan, std::map<size_t, double>>{};
  size_t n = 10;
  size_t dim = 10;
  std::vector<std::vector<double>> vx(n, std::vector<double>(dim));
  std::vector<std::map<size_t, double>> my(n);
  std::vector<DistrArraySpan> cx;
  cx.reserve(n);
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm_global(), &mpi_rank);
  MPI_Comm_size(comm_global(), &mpi_size);
  for (size_t i = 0; i < n; i++) {
    std::iota(vx[i].begin(), vx[i].end(), i + 0.5);
    my[i] = std::map<size_t, double>{{1, i + 1.0}, {3, i + 2.0}, {6, i + 3.0}, {9, i + 4.0}};
    auto crange = make_distribution_spread_remainder<size_t>(dim, mpi_size).range(mpi_rank);
    auto clength = crange.second - crange.first;
    cx.emplace_back(dim, Span<DistrArraySpan::value_type>(&vx[i][crange.first], clength), comm_global());
  }
  std::vector<double> vref(n * n), vgemm(n * n);
  std::pair<size_t, size_t> mat_dim = std::make_pair(n, n);
  Matrix<double> gemm_dot(vref, mat_dim);
  Matrix<double> ref_dot(vgemm, mat_dim);
  gemm_dot = handler.gemm_inner(cwrap(cx), cwrap(my));
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      ref_dot(i, j) = handler.dot(cx[i], my[j]);
    }
  }
  EXPECT_THAT(vgemm, Pointwise(DoubleEq(), vref));
}

TEST(TestGemm, ddisksparse_inner) {
  using molpro::mpi::comm_global;

  auto handler = ArrayHandlerDDiskSparse<DistrArrayFile, std::map<size_t, double>>{};
  size_t n = 10;
  size_t dim = 10;
  std::vector<std::vector<double>> vx(n, std::vector<double>(dim));
  std::vector<std::map<size_t, double>> my(n);
  std::vector<DistrArrayFile> cx;
  cx.reserve(n);
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm_global(), &mpi_rank);
  MPI_Comm_size(comm_global(), &mpi_size);
  for (size_t i = 0; i < n; i++) {
    std::iota(vx[i].begin(), vx[i].end(), i + 0.5);
    my[i] = std::map<size_t, double>{{1, i + 1.0}, {3, i + 2.0}, {6, i + 3.0}, {9, i + 4.0}};
    cx.emplace_back(dim);
    auto crange = make_distribution_spread_remainder<size_t>(dim, mpi_size).range(mpi_rank);
    cx.back().put(crange.first, crange.second, &(*(vx[i].cbegin() + crange.first)));
  }
  std::vector<double> vref(n * n), vgemm(n * n);
  std::pair<size_t, size_t> mat_dim = std::make_pair(n, n);
  Matrix<double> gemm_dot(vref, mat_dim);
  Matrix<double> ref_dot(vgemm, mat_dim);
  gemm_dot = handler.gemm_inner(cwrap(cx), cwrap(my));
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      ref_dot(i, j) = handler.dot(cx[i], my[j]);
    }
  }
  EXPECT_THAT(vgemm, Pointwise(DoubleEq(), vref));
}

TEST(TestGemm, distr_outer) {
  using molpro::mpi::comm_global;

  auto handler = ArrayHandlerDistr<DistrArraySpan, DistrArraySpan>{};
  size_t n = 10;
  size_t dim = 10;
  std::vector<std::vector<double>> vx(n, std::vector<double>(dim)), vy(n, std::vector<double>(dim)),
      vz(n, std::vector<double>(dim));
  std::vector<DistrArraySpan> cx, cy, cz;
  cx.reserve(n);
  cy.reserve(n);
  cz.reserve(n);
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm_global(), &mpi_rank);
  MPI_Comm_size(comm_global(), &mpi_size);
  for (size_t i = 0; i < n; i++) {
    std::iota(vx[i].begin(), vx[i].end(), i + 0.5);
    std::iota(vy[i].begin(), vy[i].end(), i + 0.5);
    std::iota(vz[i].begin(), vz[i].end(), i + 0.5);
    auto crange = make_distribution_spread_remainder<size_t>(dim, mpi_size).range(mpi_rank);
    auto clength = crange.second - crange.first;
    cx.emplace_back(dim, Span<DistrArraySpan::value_type>(&vx[i][crange.first], clength), comm_global());
    cy.emplace_back(dim, Span<DistrArraySpan::value_type>(&vy[i][crange.first], clength), comm_global());
    cz.emplace_back(dim, Span<DistrArraySpan::value_type>(&vz[i][crange.first], clength), comm_global());
  }
  std::vector<double> coeff(n * n);
  std::iota(coeff.begin(), coeff.end(), 1);
  std::pair<size_t, size_t> mat_dim = std::make_pair(n, n);
  Matrix<double> alpha(coeff, mat_dim);
  handler.gemm_outer(alpha, cwrap(cx), wrap(cy));
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      handler.axpy(alpha(i, j), cx[i], cz[j]);
    }
  }
  for (size_t i = 0; i < n; i++) {
    EXPECT_THAT(vy[i], Pointwise(DoubleEq(), vz[i]));
  }
}

TEST(TestGemm, distrddisk_outer) {
  using molpro::mpi::comm_global;

  auto handler = ArrayHandlerDistrDDisk<DistrArraySpan, DistrArrayFile>{};
  size_t n = 250;
  size_t dim = 250;
  std::vector<std::vector<double>> vx(n, std::vector<double>(dim)), vy(n, std::vector<double>(dim)),
      vz(n, std::vector<double>(dim));
  std::vector<DistrArraySpan> cx, cy;
  std::vector<DistrArrayFile> cz;
  cx.reserve(n);
  cy.reserve(n);
  cz.reserve(n);
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm_global(), &mpi_rank);
  MPI_Comm_size(comm_global(), &mpi_size);
  for (size_t i = 0; i < n; i++) {
    std::iota(vx[i].begin(), vx[i].end(), i + 0.5);
    std::iota(vy[i].begin(), vy[i].end(), i + 0.5);
    std::iota(vz[i].begin(), vz[i].end(), i + 0.5);
    cz.emplace_back(dim);
    auto crange = make_distribution_spread_remainder<size_t>(dim, mpi_size).range(mpi_rank);
    auto clength = crange.second - crange.first;
    cx.emplace_back(dim, Span<DistrArraySpan::value_type>(&vx[i][crange.first], clength), comm_global());
    cy.emplace_back(dim, Span<DistrArraySpan::value_type>(&vy[i][crange.first], clength), comm_global());
    cz.back().put(crange.first, crange.second, &(*(vz[i].cbegin() + crange.first)));
  }
  for (size_t i = 0; i != cz.size(); i++) {
    cz[i].set_buffer_size(64);
  }
  std::vector<double> coeff(n * n);
  std::iota(coeff.begin(), coeff.end(), 1);
  std::pair<size_t, size_t> mat_dim = std::make_pair(n, n);
  Matrix<double> alpha(coeff, mat_dim);
  handler.gemm_outer(alpha, cwrap(cz), wrap(cx));
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      handler.axpy(alpha(i, j), cz[i], cy[j]);
    }
  }
  for (size_t i = 0; i < n; i++) {
    EXPECT_THAT(vy[i], Pointwise(DoubleEq(), vx[i]));
  }
}

TEST(TestGemm, ddiskdistr_outer) {
  using molpro::mpi::comm_global;

  auto handler = ArrayHandlerDDiskDistr<DistrArrayFile, DistrArraySpan>{};
  size_t n = 10;
  size_t dim = 10;
  std::vector<std::vector<double>> vx(n, std::vector<double>(dim)), vy(n, std::vector<double>(dim)),
      vz(n, std::vector<double>(dim));
  std::vector<DistrArrayFile> cx, cy;
  std::vector<DistrArraySpan> cz;
  cx.reserve(n);
  cy.reserve(n);
  cz.reserve(n);
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm_global(), &mpi_rank);
  MPI_Comm_size(comm_global(), &mpi_size);
  for (size_t i = 0; i < n; i++) {
    std::iota(vx[i].begin(), vx[i].end(), i + 0.5);
    std::iota(vy[i].begin(), vy[i].end(), i + 0.5);
    std::iota(vz[i].begin(), vz[i].end(), i + 0.5);
    cx.emplace_back(dim);
    cy.emplace_back(dim);
    auto crange = make_distribution_spread_remainder<size_t>(dim, mpi_size).range(mpi_rank);
    auto clength = crange.second - crange.first;
    cx.back().put(crange.first, crange.second, &(*(vx[i].cbegin() + crange.first)));
    cy.back().put(crange.first, crange.second, &(*(vy[i].cbegin() + crange.first)));
    cz.emplace_back(dim, Span<DistrArraySpan::value_type>(&vz[i][crange.first], clength), comm_global());
  }
  std::vector<double> coeff(n * n);
  std::iota(coeff.begin(), coeff.end(), 1);
  std::pair<size_t, size_t> mat_dim = std::make_pair(n, n);
  Matrix<double> alpha(coeff, mat_dim);
  // the arrays cz in memory are addressed in place, while cx on disk are updated page by page
  cx.front().set_buffer_size(3);
  handler.gemm_outer(alpha, cwrap(cz), wrap(cx));
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      handler.axpy(alpha(i, j), cz[i], cy[j]);
    }
  }
  for (size_t i = 0; i < n; i++) {
    std::vector<double> tx(dim, 0), ty(dim, 0);
    auto crange = make_distribution_spread_remainder<size_t>(dim, mpi_size).range(mpi_rank);
    cx[i].get(crange.first, crange.second, &(*(tx.begin() + crange.first)));
    cy[i].get(crange.first, crange.second, &(*(ty.begin() + crange.first)));
    EXPECT_THAT(ty, Pointwise(DoubleEq(), tx));
  }
}

TEST(TestGemm, ddisk_outer) {
  using molpro::mpi::comm_global;

  auto handler = ArrayHandlerDDisk<DistrArrayFile, DistrArrayFile>{};
  size_t n = 10;
  size_t dim = 10;
  std::vector<std::vector<double>> vx(n, std::vector<double>(dim)), vy(n, std::vector<double>(dim)),
      vz(n, std::vector<double>(dim));
  std::vector<DistrArrayFile> cx, cy, cz;
  cx.reserve(n);
  cy.reserve(n);
  cz.reserve(n);
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm_global(), &mpi_rank);
  MPI_Comm_size(comm_global(), &mpi_size);
  for (size_t i = 0; i < n; i++) {
    std::iota(vx[i].begin(), vx[i].end(), i + 0.5);
    std::iota(vy[i].begin(), vy[i].end(), i + 0.5);
    std::iota(vz[i].begin(), vz[i].end(), i + 0.5);
    cx.emplace_back(dim);
    cy.emplace_back(dim);
    cz.emplace_back(dim);
    auto crange = make_distribution_spread_remainder<size_t>(dim, mpi_size).range(mpi_rank);
    cx.back().put(crange.first, crange.second, &(*(vx[i].cbegin() + crange.first)));
    cy.back().put(crange.first, crange.second, &(*(vy[i].cbegin() + crange.first)));
    cz.back().put(crange.first, crange.second, &(*(vz[i].cbegin() + crange.first)));
  }
  std::vector<double> coeff(n * n);
  std::iota(coeff.begin(), coeff.end(), 1);
  std::pair<size_t, size_t> mat_dim = std::make_pair(n, n);
  Matrix<double> alpha(coeff, mat_dim);
  handler.gemm_outer(alpha, cwrap(cz), wrap(cx));
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      handler.axpy(alpha(i, j), cz[i], cy[j]);
    }
  }
  for (size_t i = 0; i < n; i++) {
    std::vector<double> tx(dim, 0), ty(dim, 0);
    auto crange = cx.back().distribution().range(mpi_rank);
    cx[i].get(crange.first, crange.second, &(*(tx.begin() + crange.first)));
    cy[i].get(crange.first, crange.second, &(*(ty.begin() + crange.first)));
    EXPECT_THAT(ty, Pointwise(DoubleEq(), tx));
  }
}

TEST(TestGemm, distrsparse_outer) {
  using molpro::mpi::comm_global;

  auto handler = ArrayHandlerDistrSparse<DistrArraySpan, std::map<size_t, double>>{};
  size_t n = 10;
  size_t dim = 10;
  std::vector<std::vector<double>> vx(n, std::vector<double>(dim)), vy(n, std::vector<double>(dim));
  std::vector<std::map<size_t, double>> my(n);
  std::vector<DistrArraySpan> cx, cy;
  cx.reserve(n);
  cy.reserve(n);
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm_global(), &mpi_rank);
  MPI_Comm_size(comm_global(), &mpi_size);
  for (size_t i = 0; i < n; i++) {
    std::iota(vx[i].begin(), vx[i].end(), i + 0.5);
    std::iota(vy[i].begin(), vy[i].end(), i + 0.5);
    my[i] = std::map<size_t, double>{{1, i + 1.0}, {3, i + 2.0}, {6, i + 3.0}, {9, i + 4.0}};
    auto crange = make_distribution_spread_remainder<size_t>(dim, mpi_size).range(mpi_rank);
    auto clength = crange.second - crange.first;
    cx.emplace_back(dim, Span<DistrArraySpan::value_type>(&vx[i][crange.first], clength), comm_global());
    cy.emplace_back(dim, Span<DistrArraySpan::value_type>(&vy[i][crange.first], clength), comm_global());
  }
  std::vector<double> coeff(n * n);
  std::iota(coeff.begin(), coeff.end(), 1);
  std::pair<size_t, size_t> mat_dim = std::make_pair(n, n);
  Matrix<double> alpha(coeff, mat_dim);
  handler.gemm_outer(alpha, cwrap(my), wrap(cx));
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      handler.axpy(alpha(j, i), my[j], cy[i]);
    }
  }
  for (size_t i = 0; i < n; i++) {
    EXPECT_THAT(vy[i], Pointwise(DoubleEq(), vx[i]));
  }
}

TEST(TestGemm, ddisksparse_outer) {
  using molpro::mpi::comm_global;

  auto handler = ArrayHandlerDDiskSparse<DistrArrayFile, std::map<size_t, double>>{};
  size_t n = 10;
  size_t dim = 10;
  std::vector<std::vector<double>> vx(n, std::vector<double>(dim)), vy(n, std::vector<double>(dim));
  std::vector<std::map<size_t, double>> my(n);
  std::vector<DistrArrayFile> cx, cy;
  cx.reserve(n);
  cy.reserve(n);
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm_global(), &mpi_rank);
  MPI_Comm_size(comm_global(), &mpi_size);
  for (size_t i = 0; i < n; i++) {
    std::iota(vx[i].begin(), vx[i].end(), i + 0.5);
    std::iota(vy[i].begin(), vy[i].end(), i + 0.5);
    my[i] = std::map<size_t, double>{{1, i + 1.0}, {3, i + 2.0}, {6, i + 3.0}, {9, i + 4.0}};
    cx.emplace_back(dim);
    cy.emplace_back(dim);
    auto crange = cx.back().distribution().range(mpi_rank);
    cx.back().put(crange.first, crange.second, &(*(vx[i].cbegin() + crange.first)));
    cy.back().put(crange.first, crange.second, &(*(vy[i].cbegin() + crange.first)));
  }
  std::vector<double> coeff(n * n);
  std::iota(coeff.begin(), coeff.end(), 1);
  std::pair<size_t, size_t> mat_dim = std::make_pair(n, n);
  Matrix<double> alpha(coeff, mat_dim);
  handler.gemm_outer(alpha, cwrap(my), wrap(cx));
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      handler.axpy(alpha(j, i), my[j], cy[i]);
    }
  }
  for (size_t i = 0; i < n; i++) {
    std::vector<double> tx(dim, 0), ty(dim, 0);
    auto crange = cx.back().distribution().range(mpi_rank);
    cx[i].get(crange.first, crange.second, &(*(tx.begin() + crange.first)));
    cy[i].get(crange.first, crange.second, &(*(ty.begin() + crange.first)));
    EXPECT_THAT(ty, Pointwise(DoubleEq(), tx));
  }
}

class BufferedDistrArrayFileTest : public testing::TestWithParam<size_t> {};

TEST_P(BufferedDistrArrayFileTest, GEMM) {
  using molpro::mpi::comm_global;

  const size_t n = GetParam();

  auto handler = ArrayHandlerDistrDDisk<DistrArraySpan, DistrArrayFile>{};
  size_t dim = 10000; // height
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm_global(), &mpi_rank);
  MPI_Comm_size(comm_global(), &mpi_size);

  std::vector<double> coeff(n * n);
  std::iota(coeff.begin(), coeff.end(), 1);
  std::pair<size_t, size_t> mat_dim = std::make_pair(n, n);
  Matrix<double> alpha(coeff, mat_dim);

  for (size_t stride_multiplier = 0; stride_multiplier < 2; ++stride_multiplier) {
    auto [cx, cy, cz] = molpro::linalg::test::get_contiguous(n, dim);
    decltype(cx) cx_selection;
    for (size_t i = 0; i < n; i += stride_multiplier * i + 1)
      cx_selection.emplace_back(cx[i]);
    decltype(alpha) alpha_selection({n, cx_selection.size()});
    {
      size_t offset = 0;
      for (size_t i = 0; i < n; i += stride_multiplier * i + 1) {
        for (size_t j = 0; j < n; ++j)
          alpha_selection(j, offset) = alpha(j, i);
        ++offset;
      }
      assert(offset == cx_selection.size());
    }

    Matrix<double> expected_result({cx_selection.size(), dim});
    const auto distribution =
        cx_selection.empty() ? molpro::linalg::array::util::make_distribution_spread_remainder<size_t>(dim, mpi_size)
                             : cx_selection.front().distribution();
    const auto range = distribution.range(mpi_rank);

    for (size_t i = 0; i < cx_selection.size(); i++) {
      for (size_t j = range.first; j < range.second; j++) {
        expected_result(i, j) = cx_selection[i][j];
      }
    }
    for (size_t i = 0; i < n; i++) {
      std::vector<double> czbuf(range.second - range.first);
      cz[i].get(range.first, range.second, czbuf.data());
      for (size_t j = range.first; j < range.second; j++) {
        for (size_t k = 0; k < cx_selection.size(); ++k) {
          expected_result(k, j) += czbuf[j - range.first] * alpha_selection(i, k);
        }
      }
    }
#ifdef HAVE_MPI_H
    for (size_t i = 0; i < cx_selection.size(); i++) {
      for (int rank = 0; rank < mpi_size; ++rank) {
        MPI_Bcast(&(expected_result(i, distribution.range(rank).first)),
                  distribution.range(rank).second - distribution.range(rank).first, MPI_DOUBLE, rank,
                  comm_global());
      }
    }
#endif

    if (cx_selection.size() != cx.size()) {
      EXPECT_THROW(handler.gemm_outer(alpha, cwrap(cz), wrap(cx_selection)), std::out_of_range);
    }

    handler.gemm_outer(alpha_selection, cwrap(cz), wrap(cx_selection));

    Matrix<double> actual_result({cx_selection.size(), dim});
    for (size_t i = 0; i < cx_selection.size(); i++) {
      for (size_t j = range.first; j < range.second; j++)
        actual_result(i, j) = cx_selection[i][j];
#ifdef HAVE_MPI_H
     for (int rank = 0; rank < mpi_size; ++rank)
       MPI_Bcast(&(actual_result(i, distribution.range(rank).first)),
                 distribution.range(rank).second - distribution.range(rank).first, MPI_DOUBLE, rank,
                 comm_global());
#endif
    }
    if (false) {
      molpro::linalg::array::util::LockMPI3 lock(comm_global());
      std::cout << "expected_result" << as_string(expected_result) << std::endl;
      std::cout << "actual result" << as_string(actual_result) << std::endl;
    }
    EXPECT_THAT(actual_result.data(), Pointwise(DoubleEq(), expected_result.data()));

    auto alpha_selection_calculated = handler.gemm_inner(cwrap(cx_selection), cwrap(cz));
    auto alpha_selection_expected = alpha_selection_calculated;
    alpha_selection_expected.fill(0);
    for (size_t i = 0; i < cx_selection.size(); ++i) {
      for (size_t k = 0; k < n; k++) {
        std::vector<double> czbuf(range.second - range.first);
        cz[k].get(range.first, range.second, czbuf.data());
        for (size_t j = range.first; j < range.second; j++) {
          //                        std::cout << i<<" "<<j<<" "<<k<<" "<<czbuf[j-range.first]<<"
          //                        "<<cx_selection[i][j]<<std::endl;
          alpha_selection_expected(i, k) += czbuf[j - range.first] * cx_selection[i][j];
        }
      }
    }
#ifdef HAVE_MPI_H
    MPI_Allreduce(MPI_IN_PLACE, const_cast<double*>(alpha_selection_expected.data().data()),
                  alpha_selection_expected.size(), MPI_DOUBLE, MPI_SUM, comm_global());
#endif
    //      std::cout<< "Calculated alpha:" << as_string(alpha_selection_calculated)
    //                             << "\nExpected alpha:" << as_string(alpha_selection_expected)<<std::endl;

    if (alpha_selection_expected.size() > 0) {
      EXPECT_THAT(
          alpha_selection_calculated.data(),
          Pointwise(::testing::DoubleNear(1e-12 * alpha_selection_expected(0, 0)), alpha_selection_expected.data()))
          << "Calculated alpha:" << as_string(alpha_selection_calculated)
          << "\nExpected alpha:" << as_string(alpha_selection_expected)
          << "\nstride_multiplier = " << stride_multiplier;
    }
  }
}

INSTANTIATE_TEST_SUITE_P(,BufferedDistrArrayFileTest, testing::Range(static_cast<size_t>(0), static_cast<size_t>(9)));

// The buffered gemm must serve every array held on disk, not only DistrArrayFile itself (issue #119). A subclass of
// DistrArrayFile stands in for another DistrArrayDisk, such as DistrArrayHDF5: the buffered path reads it in chunks
// through get(), whereas the unbuffered path would read whole vectors through local_buffer().
namespace {
struct CountingFile : DistrArrayFile {
  using DistrArrayFile::DistrArrayFile;
  static inline int local_buffer_calls = 0;
  [[nodiscard]] std::unique_ptr<LocalBuffer> local_buffer() override {
    ++local_buffer_calls;
    return DistrArrayFile::local_buffer();
  }
  [[nodiscard]] std::unique_ptr<const LocalBuffer> local_buffer() const override {
    ++local_buffer_calls;
    return DistrArrayFile::local_buffer();
  }
};
} // namespace

TEST(TestGemm, distrarraydisk_subclass_buffered) {
  using molpro::linalg::array::util::gemm_inner_distr_distr;
  using molpro::linalg::array::util::gemm_outer_distr_distr;
  using molpro::mpi::comm_global;
  static_assert(molpro::linalg::array::util::is_disk_array_v<CountingFile>);
  const size_t n = 4, dim = 30;
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm_global(), &mpi_rank);
  MPI_Comm_size(comm_global(), &mpi_size);
  const auto crange = make_distribution_spread_remainder<size_t>(dim, mpi_size).range(mpi_rank);
  const auto clength = crange.second - crange.first;
  std::vector<std::vector<double>> vx(n, std::vector<double>(dim)), vy(n, std::vector<double>(dim)),
      vz(n, std::vector<double>(dim)), vzref(n, std::vector<double>(dim));
  std::vector<CountingFile> cx;
  std::vector<DistrArrayFile> cxref;
  std::vector<DistrArraySpan> cy, cz, czref;
  cx.reserve(n);
  cxref.reserve(n);
  cy.reserve(n);
  cz.reserve(n);
  czref.reserve(n);
  for (size_t i = 0; i < n; i++) {
    for (size_t k = 0; k < dim; ++k) {
      vx[i][k] = std::cos(1.0 + i + 0.3 * k);
      vy[i][k] = std::sin(2.0 + i + 0.7 * k);
      vz[i][k] = vzref[i][k] = 0.1 * (i + k);
    }
    cx.emplace_back(dim);
    cx.back().put(crange.first, crange.second, &vx[i][crange.first]);
    cxref.emplace_back(dim);
    cxref.back().put(crange.first, crange.second, &vx[i][crange.first]);
    cy.emplace_back(dim, Span<double>(&vy[i][crange.first], clength), comm_global());
    cz.emplace_back(dim, Span<double>(&vz[i][crange.first], clength), comm_global());
    czref.emplace_back(dim, Span<double>(&vzref[i][crange.first], clength), comm_global());
  }
  std::vector<double> coeff(n * n);
  std::iota(coeff.begin(), coeff.end(), 1);
  const Matrix<double> alpha(coeff, std::make_pair(n, n));

  CountingFile::local_buffer_calls = 0;
  const auto inner_yx = gemm_inner_distr_distr(cwrap(cy), cwrap(cx));
  const auto inner_xy = gemm_inner_distr_distr(cwrap(cx), cwrap(cy));
  gemm_outer_distr_distr(alpha, cwrap(cx), wrap(cz));
  EXPECT_EQ(CountingFile::local_buffer_calls, 0);

  // the same operations on DistrArrayFile, and computed directly from the data
  const auto ref_yx = gemm_inner_distr_distr(cwrap(cy), cwrap(cxref));
  const auto ref_xy = gemm_inner_distr_distr(cwrap(cxref), cwrap(cy));
  gemm_outer_distr_distr(alpha, cwrap(cxref), wrap(czref));
  Matrix<double> direct_yx({n, n}), direct_xy({n, n});
  std::vector<std::vector<double>> direct_z(n, std::vector<double>(dim));
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j) {
      direct_yx(i, j) = std::inner_product(vy[i].begin(), vy[i].end(), vx[j].begin(), 0.0);
      direct_xy(i, j) = std::inner_product(vx[i].begin(), vx[i].end(), vy[j].begin(), 0.0);
    }
  for (size_t j = 0; j < n; ++j)
    for (size_t k = 0; k < dim; ++k) {
      direct_z[j][k] = 0.1 * (j + k);
      for (size_t i = 0; i < n; ++i)
        direct_z[j][k] += alpha(i, j) * vx[i][k];
    }
  for (size_t i = 0; i < n * n; ++i) {
    EXPECT_NEAR(inner_yx.data()[i], direct_yx.data()[i], 1e-12) << "inner (in memory, disk) " << i;
    EXPECT_NEAR(inner_xy.data()[i], direct_xy.data()[i], 1e-12) << "inner (disk, in memory) " << i;
    EXPECT_NEAR(ref_yx.data()[i], direct_yx.data()[i], 1e-12) << "DistrArrayFile inner (in memory, disk) " << i;
    EXPECT_NEAR(ref_xy.data()[i], direct_xy.data()[i], 1e-12) << "DistrArrayFile inner (disk, in memory) " << i;
  }
  for (size_t j = 0; j < n; ++j)
    for (size_t k = crange.first; k < crange.second; ++k) {
      EXPECT_NEAR(vz[j][k], direct_z[j][k], 1e-12) << "outer " << j << ", " << k;
      EXPECT_NEAR(vzref[j][k], direct_z[j][k], 1e-12) << "DistrArrayFile outer " << j << ", " << k;
    }
}
