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
#ifdef LINEARALGEBRA_ARRAY_MPI3
#include <molpro/linalg/array/DistrArrayMPI3.h>
#endif
#include <molpro/linalg/array/util.h>
#include <molpro/linalg/array/util/gemm.h>
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
#ifdef LINEARALGEBRA_ARRAY_MPI3
using molpro::linalg::array::DistrArrayMPI3;
#endif
using molpro::linalg::array::Span;
using molpro::linalg::array::util::LockMPI3;
using molpro::linalg::array::util::ScopeLock;
//using molpro::linalg::test::mpi_comm;
using molpro::linalg::itsolv::wrap;
using molpro::linalg::itsolv::cwrap;
using molpro::linalg::itsolv::subspace::Matrix;

using ::testing::ContainerEq;
using ::testing::DoubleEq;
using ::testing::Each;
using ::testing::Pointwise;

TEST(TestGemmInner, distrdistr) {
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
    //cx.emplace_back(dim, mpi_comm);
    //cy.emplace_back(dim, mpi_comm);
    cx.emplace_back(dim);
    cy.emplace_back(dim);
    auto crange = cx.back().distribution().range(mpi_rank);
    auto clength = crange.second - crange.first;
    cx.back().allocate_buffer(Span<DistrArraySpan::value_type>(&vx[i][crange.first], clength));
    cy.back().allocate_buffer(Span<DistrArraySpan::value_type>(&vx[i][crange.first], clength));
  }
  std::vector<double> vref(n*n), vgemm(n*n);
  //std::iota(coeff.begin(), coeff.end(), 1);
  std::pair<size_t,size_t> mat_dim = std::make_pair(n,n);
  //Matrix<double> alpha(coeff, mat_dim);
  Matrix<double> gemm_dot(vref, mat_dim);
  Matrix<double> ref_dot(vgemm, mat_dim);
  gemm_dot = handler.gemm_inner(cwrap(cx),cwrap(cy));
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      ref_dot(i, j) = handler.dot(cx[i], cy[j]);
    }
  }
  std::cout << "ref dot:" << std::endl;
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      ref_dot(i, j) = handler.dot(cx[i], cy[j]);
      std::cout << ref_dot(i, j) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "gemm dot:" << std::endl;
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      std::cout << gemm_dot(i, j) << " ";
    }
    std::cout << std::endl;
  }
  EXPECT_THAT(vgemm, Pointwise(DoubleEq(), vref));
}
