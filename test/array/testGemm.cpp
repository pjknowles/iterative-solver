#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <vector>

#include "parallel_util.h"

#include <molpro/linalg/array/ArrayHandlerDistr.h>
#include <molpro/linalg/array/ArrayHandlerDDisk.h>
#include <molpro/linalg/array/ArrayHandlerSparse.h>
#include <molpro/linalg/array/ArrayHandlerIterable.h>
#include <molpro/linalg/array/ArrayHandlerDistrDDisk.h>
#include <molpro/linalg/array/ArrayHandlerDDiskDistr.h>
#include <molpro/linalg/array/ArrayHandlerDistrSparse.h>
#include <molpro/linalg/array/ArrayHandlerDDiskSparse.h>
#include <molpro/linalg/array/ArrayHandlerIterableSparse.h>
#include <molpro/linalg/array/DistrArraySpan.h>
#include <molpro/linalg/array/DistrArrayFile.h>
#include <molpro/linalg/array/DistrArrayDisk.h>
#ifdef LINEARALGEBRA_ARRAY_MPI3
#include <molpro/linalg/array/DistrArrayMPI3.h>
#endif
#include <molpro/linalg/array/util.h>
#include <molpro/linalg/itsolv/wrap.h>
#include <molpro/linalg/array/util/Distribution.h>
#include <molpro/linalg/array/Span.h>
#include <molpro/linalg/itsolv/subspace/Matrix.h>

using molpro::linalg::array::ArrayHandlerDistr;
using molpro::linalg::array::ArrayHandlerDDisk;
using molpro::linalg::array::ArrayHandlerSparse;
using molpro::linalg::array::ArrayHandlerIterable;
using molpro::linalg::array::ArrayHandlerDistrDDisk;
using molpro::linalg::array::ArrayHandlerDDiskDistr;
using molpro::linalg::array::ArrayHandlerDistrSparse;
using molpro::linalg::array::ArrayHandlerDDiskSparse;
using molpro::linalg::array::ArrayHandlerIterableSparse;
using molpro::linalg::array::DistrArraySpan;
using molpro::linalg::array::DistrArrayFile;
using molpro::linalg::array::DistrArrayDisk;
#ifdef LINEARALGEBRA_ARRAY_MPI3
using molpro::linalg::array::DistrArrayMPI3;
#endif
using molpro::linalg::array::Span;
using molpro::linalg::array::util::LockMPI3;
using molpro::linalg::array::util::ScopeLock;
using molpro::linalg::itsolv::wrap;
using molpro::linalg::itsolv::cwrap;
using molpro::linalg::itsolv::subspace::Matrix;
using molpro::linalg::array::util::make_distribution_spread_remainder;

using ::testing::ContainerEq;
using ::testing::DoubleEq;
using ::testing::Each;
using ::testing::Pointwise;

TEST(TestGemm, distr_inner) {
  auto handler = ArrayHandlerDistr<DistrArraySpan,DistrArraySpan>{};
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
  std::vector<double> vref(n*n), vgemm(n*n);
  std::pair<size_t,size_t> mat_dim = std::make_pair(n,n);
  Matrix<double> gemm_dot(vref, mat_dim);
  Matrix<double> ref_dot(vgemm, mat_dim);
  gemm_dot = handler.gemm_inner(cwrap(cx),cwrap(cy));
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      ref_dot(i, j) = handler.dot(cx[i], cy[j]);
    }
  }
  EXPECT_THAT(vgemm, Pointwise(DoubleEq(), vref));
}

TEST(TestGemm, distrarrayfile_inner) {
  auto handler = ArrayHandlerDistr<DistrArraySpan,DistrArrayFile>{};
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
    cy.back().put(crange.first, crange.second ,vy[i].data());
  }
  std::vector<double> vref, vgemm;
  Matrix<double> gemm_dot({n,n});
  gemm_dot = handler.gemm_inner(cwrap(cx),cwrap(cy));
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      vgemm.push_back(gemm_dot(i,j));
      vref.push_back(handler.dot(cx[i], cy[j]));
    }
  }
  EXPECT_THAT(vgemm, Pointwise(DoubleEq(), vref));
}

TEST(TestGemm, ddiskdistr_inner) {
  auto handler = ArrayHandlerDDiskDistr<DistrArrayFile,DistrArraySpan>{};
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
  std::vector<double> vref(n*n), vgemm(n*n);
  std::pair<size_t,size_t> mat_dim = std::make_pair(n,n);
  Matrix<double> gemm_dot(vref, mat_dim);
  Matrix<double> ref_dot(vgemm, mat_dim);
  gemm_dot = handler.gemm_inner(cwrap(cx),cwrap(cy));
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      ref_dot(i, j) = handler.dot(cx[i], cy[j]);
    }
  }
  EXPECT_THAT(vgemm, Pointwise(DoubleEq(), vref));
}

TEST(TestGemm, distrddisk_inner) {
  auto handler = ArrayHandlerDistrDDisk<DistrArraySpan,DistrArrayFile>{};
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
  std::vector<double> vref(n*n), vgemm(n*n);
  std::pair<size_t,size_t> mat_dim = std::make_pair(n,n);
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
  auto handler = ArrayHandlerDDisk<DistrArrayFile,DistrArrayFile>{};
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
  std::vector<double> vref(n*n), vgemm(n*n);
  std::pair<size_t,size_t> mat_dim = std::make_pair(n,n);
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
  auto handler = ArrayHandlerDistrSparse<DistrArraySpan,std::map<size_t, double>>{};
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
    my[i] = std::map<size_t, double>{{1, i+1.0}, {3, i+2.0}, {6, i+3.0}, {9, i+4.0}};
    auto crange = make_distribution_spread_remainder<size_t>(dim, mpi_size).range(mpi_rank);
    auto clength = crange.second - crange.first;
    cx.emplace_back(dim, Span<DistrArraySpan::value_type>(&vx[i][crange.first], clength), comm_global());
  }
  std::vector<double> vref(n*n), vgemm(n*n);
  std::pair<size_t,size_t> mat_dim = std::make_pair(n,n);
  Matrix<double> gemm_dot(vref, mat_dim);
  Matrix<double> ref_dot(vgemm, mat_dim);
  gemm_dot = handler.gemm_inner(cwrap(cx),cwrap(my));
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      ref_dot(i, j) = handler.dot(cx[i], my[j]);
    }
  }
  EXPECT_THAT(vgemm, Pointwise(DoubleEq(), vref));
}

TEST(TestGemm, ddisksparse_inner) {
  auto handler = ArrayHandlerDDiskSparse<DistrArrayFile,std::map<size_t, double>>{};
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
    my[i] = std::map<size_t, double>{{1, i+1.0}, {3, i+2.0}, {6, i+3.0}, {9, i+4.0}};
    cx.emplace_back(dim);
    auto crange = make_distribution_spread_remainder<size_t>(dim, mpi_size).range(mpi_rank);
    cx.back().put(crange.first, crange.second, &(*(vx[i].cbegin() + crange.first)));
  }
  std::vector<double> vref(n*n), vgemm(n*n);
  std::pair<size_t,size_t> mat_dim = std::make_pair(n,n);
  Matrix<double> gemm_dot(vref, mat_dim);
  Matrix<double> ref_dot(vgemm, mat_dim);
  gemm_dot = handler.gemm_inner(cwrap(cx),cwrap(my));
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      ref_dot(i, j) = handler.dot(cx[i], my[j]);
    }
  }
  EXPECT_THAT(vgemm, Pointwise(DoubleEq(), vref));
}

TEST(TestGemm, distr_outer) {
  auto handler = ArrayHandlerDistr<DistrArraySpan,DistrArraySpan>{};
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
  std::vector<double> coeff(n*n);
  std::iota(coeff.begin(), coeff.end(), 1);
  std::pair<size_t,size_t> mat_dim = std::make_pair(n,n);
  Matrix<double> alpha(coeff, mat_dim);
  handler.gemm_outer(alpha, cwrap(cx),wrap(cy));
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
  auto handler = ArrayHandlerDistrDDisk<DistrArraySpan,DistrArrayFile>{};
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
  for(size_t i = 0; i != cz.size(); i++){
    cz[i].set_buffer_size(64);
  }
  std::vector<double> coeff(n*n);
  std::iota(coeff.begin(), coeff.end(), 1);
  std::pair<size_t,size_t> mat_dim = std::make_pair(n,n);
  Matrix<double> alpha(coeff, mat_dim);
  handler.gemm_outer(alpha, cwrap(cz),wrap(cx));
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
  auto handler = ArrayHandlerDDiskDistr<DistrArrayFile,DistrArraySpan>{};
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
  std::vector<double> coeff(n*n);
  std::iota(coeff.begin(), coeff.end(), 1);
  std::pair<size_t,size_t> mat_dim = std::make_pair(n,n);
  Matrix<double> alpha(coeff, mat_dim);
  handler.gemm_outer(alpha, cwrap(cz),wrap(cx));
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
  auto handler = ArrayHandlerDDisk<DistrArrayFile,DistrArrayFile>{};
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
  std::vector<double> coeff(n*n);
  std::iota(coeff.begin(), coeff.end(), 1);
  std::pair<size_t,size_t> mat_dim = std::make_pair(n,n);
  Matrix<double> alpha(coeff, mat_dim);
  handler.gemm_outer(alpha, cwrap(cz),wrap(cx));
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
  auto handler = ArrayHandlerDistrSparse<DistrArraySpan,std::map<size_t, double>>{};
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
    my[i] = std::map<size_t, double>{{1, i+1.0}, {3, i+2.0}, {6, i+3.0}, {9, i+4.0}};
    auto crange = make_distribution_spread_remainder<size_t>(dim, mpi_size).range(mpi_rank);
    auto clength = crange.second - crange.first;
    cx.emplace_back(dim, Span<DistrArraySpan::value_type>(&vx[i][crange.first], clength), comm_global());
    cy.emplace_back(dim, Span<DistrArraySpan::value_type>(&vy[i][crange.first], clength), comm_global());
  }
  std::vector<double> coeff(n*n);
  std::iota(coeff.begin(), coeff.end(), 1);
  std::pair<size_t,size_t> mat_dim = std::make_pair(n,n);
  Matrix<double> alpha(coeff, mat_dim);
  handler.gemm_outer(alpha, cwrap(my),wrap(cx));
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
  auto handler = ArrayHandlerDDiskSparse<DistrArrayFile,std::map<size_t, double>>{};
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
    my[i] = std::map<size_t, double>{{1, i+1.0}, {3, i+2.0}, {6, i+3.0}, {9, i+4.0}};
    cx.emplace_back(dim);
    cy.emplace_back(dim);
    auto crange = cx.back().distribution().range(mpi_rank);
    cx.back().put(crange.first, crange.second, &(*(vx[i].cbegin() + crange.first)));
    cy.back().put(crange.first, crange.second, &(*(vy[i].cbegin() + crange.first)));
  }
  std::vector<double> coeff(n*n);
  std::iota(coeff.begin(), coeff.end(), 1);
  std::pair<size_t,size_t> mat_dim = std::make_pair(n,n);
  Matrix<double> alpha(coeff, mat_dim);
  handler.gemm_outer(alpha, cwrap(my),wrap(cx));
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

TEST(TestGemm, buffered_DistrArrayFile) {
  auto handler = ArrayHandlerDistrDDisk<DistrArraySpan, DistrArrayFile>{};
  for (size_t n = 0; n < 9; ++n) {
    size_t dim = 4099; // height
    int mpi_rank, mpi_size;
    MPI_Comm_rank(comm_global(), &mpi_rank);
    MPI_Comm_size(comm_global(), &mpi_size);

    auto [cx, cy, cz] = molpro::linalg::test::get_contiguous(n, dim);

    std::vector<double> coeff(n * n);
    std::iota(coeff.begin(), coeff.end(), 1);
    std::pair<size_t, size_t> mat_dim = std::make_pair(n, n);
    Matrix<double> alpha(coeff, mat_dim);

    for (size_t stride_multiplier = 0; stride_multiplier < 2; ++stride_multiplier) {
      decltype(cx) cx_selection;
      for (size_t i = 0; i < n; i += stride_multiplier * i + 1)
        cx_selection.emplace_back(cx[i]);
      decltype(alpha) alpha_selection({n, cx_selection.size()});
      size_t offset = 0;
      for (size_t i = 0; i < n; i += stride_multiplier * i + 1) {
        for (size_t j = 0; j < n; ++j)
          alpha_selection(j, offset) = alpha(j, i);
        ++offset;
      }
      //      std::cout << "alpha" << as_string(alpha) << std::endl;
      //      std::cout << "alpha_selection" << as_string(alpha_selection) << std::endl;

      Matrix<double> expected_result({cx_selection.size(), dim});
      for (size_t i = 0; i < cx_selection.size(); i++) {
        //        std::cout << "cx_selection";
        for (size_t j = 0; j < dim; j++) {
          expected_result(i, j) = cx_selection[i][j];
          //          std::cout << " " << cx_selection[i][j];
        }
        //        std::cout << std::endl;
      }
      //      std::cout << "expected_result before z contributions" << as_string(expected_result) << std::endl;
      for (size_t i = 0; i < n; i++) {
        //        std::cout << "cz";
        for (size_t j = 0; j < dim; j++) {
          //          std::cout << " " << cz[i][j];
          for (size_t k = 0; k < cx_selection.size(); ++k)
            expected_result(k, j) += cz[i].local_buffer()->data()[j] * alpha(i, k);
        }
        //        std::cout << std::endl;
      }

      if (cx_selection.size() != cx.size()) {
        EXPECT_THROW(handler.gemm_outer(alpha, cwrap(cz), wrap(cx_selection)), std::out_of_range);
      }
      handler.gemm_outer(alpha_selection, cwrap(cz), wrap(cx_selection));
      Matrix<double> actual_result({cx_selection.size(), dim});
      for (size_t i = 0; i < cx_selection.size(); i++)
        for (size_t j = 0; j < dim; j++)
          actual_result(i, j) = cx_selection[i][j];
      //        std::cout << "expected_result" << as_string(expected_result) << std::endl;
      //        std::cout << "actual result" << as_string(actual_result) << std::endl;
      EXPECT_THAT(actual_result.data(), Pointwise(DoubleEq(), expected_result.data()));
    }
  }
}