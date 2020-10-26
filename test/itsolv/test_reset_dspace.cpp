#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/itsolv/reset_dspace.h>

using molpro::linalg::itsolv::detail::construct_overlap_with_solutions;
using molpro::linalg::itsolv::detail::transform_PQSol_to_PQD_subspace;
using molpro::linalg::itsolv::subspace::Matrix;
using molpro::linalg::itsolv::subspace::xspace::Dimensions;
using ::testing::DoubleEq;
using ::testing::Pointwise;

TEST(reset_dspace, construct_overlap_with_solutions_empty) {
  using value_type = double;
  const auto dims = Dimensions{0, 0, 0};
  const auto ov = Matrix<value_type>{};
  const auto solutions = Matrix<value_type>{};
  auto result = construct_overlap_with_solutions(solutions, ov, dims);
  ASSERT_TRUE(result.empty());
}

TEST(reset_dspace, construct_overlap_with_solutions) {
  using value_type = double;
  const auto nP = 2, nQ = 3, nD = 1, nSol = 3;
  const auto dims = Dimensions{nP, nQ, nD};
  //@formatter:off
  const auto ov = Matrix<value_type>{{1, 0, 0, 1, 0, 1,
                                      0, 1, 0, 0, 0, 1,
                                      0, 0, 1, 0, 1, 1,
                                      1, 0, 0, 3, 1, 2,
                                      0, 0, 1, 1, 3, 2,
                                      1, 1, 1, 2, 2, 5}, {dims.nX, dims.nX}};
  const auto solutions = Matrix<value_type>{{1, 0, 0, 0, 0, 0,
                                             1, 1, 1, 1, 1, 1,
                                             0, 0, 0, 0, 0, 0}, {nSol, dims.nX}};
  const auto reference = Matrix<value_type>{{1, 0, 0, 1, 0, 1, 3, 0,
                                             0, 1, 0, 0, 0, 0, 2, 0,
                                             0, 0, 1, 0, 1, 0, 3, 0,
                                             1, 0, 0, 3, 1, 1, 7, 0,
                                             0, 0, 1, 1, 3, 0, 7, 0,
                                             1, 0, 0, 1, 0, 1, 3, 0,
                                             3, 2, 3, 7, 7, 3, 34,0,
                                             0, 0, 0, 0, 0, 0, 0, 0},
                                            {dims.nP + dims.nQ + nSol, dims.nP + dims.nQ + nSol}};
  //@formatter:on
  auto result = construct_overlap_with_solutions(solutions, ov, dims);
  ASSERT_FALSE(result.empty());
  ASSERT_EQ(result.rows(), reference.rows());
  ASSERT_EQ(result.cols(), reference.cols());
  ASSERT_THAT(result.data(), Pointwise(DoubleEq(), reference.data()));
}

TEST(reset_dspace, transform_PQSol_to_PQD_subspace){
  const auto nP = 2, nQ = 3, nD = 1, nSol = 3;
  const auto dims = Dimensions{nP, nQ, nD};
  //@formatter:off
  const auto solutions = Matrix<double>{{1, 0, 0, 0, 0, 0,
                                         1, 1, 1, 1, 1, 1,
                                         0, 0, 0, 0, 0, 0},
                                        {nSol, dims.nX}};
  const auto lin_trans_PQSol = Matrix<double>{{0, 1, 0, 1, 0, 1, 2, 3,
                                               1, 0, 1, 0, 1, 1, 1, 1,
                                               0, 0, 0, 0, 0, 0, 0, 0},
                                              {nSol, nP + nQ + nSol}};
  const auto reference_PQD = Matrix<double>{{1, 2, 1, 2, 1, 1,
                                             2, 0, 1, 0, 1, 0,
                                             0, 0, 0, 0, 0, 0},
                                            {nSol, dims.nX}};
  //@formatter:on
  const auto result = transform_PQSol_to_PQD_subspace(lin_trans_PQSol, solutions, dims);
  ASSERT_FALSE(result.empty());
  ASSERT_EQ(result.rows(), reference_PQD.rows());
  ASSERT_EQ(result.cols(), reference_PQD.cols());
  ASSERT_THAT(result.data(), Pointwise(DoubleEq(), reference_PQD.data()));
}
