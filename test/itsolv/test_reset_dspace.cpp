#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/itsolv/reset_dspace.h>

using molpro::linalg::itsolv::detail::construct_projected_solutions;
using molpro::linalg::itsolv::detail::select_full_solutions;
using molpro::linalg::itsolv::subspace::Matrix;
using molpro::linalg::itsolv::subspace::xspace::Dimensions;
using ::testing::DoubleEq;
using ::testing::Pointwise;


TEST(itsolv_reset_dspace, construct_projected_solutions) {
  const size_t nC = 2;
  const auto dims = Dimensions(1, 1, 2);
  const auto solutions = Matrix<double>({1.0, 1.0, 1.0, 1.0, //
                                         0.1, 0.2, 1.0, -1.0},
                                        {nC, dims.nX});
  const auto overlap = Matrix<double>({1, 0, 1, 0, //
                                       0, 1, 1, 0, //
                                       1, 1, 4, 1, //
                                       0, 0, 1, 1},
                                      {dims.nX, dims.nX});
  const auto ref_overlap_CC = Matrix<double>({7., 3., //
                                              3., 3.},
                                             {nC, nC});
  const auto ref_proj_solutions = Matrix<double>({1.0, 1.0, //
                                                  1.0, -1.0},
                                                 {nC, dims.nD});
  Matrix<double> res_proj_solutions, res_overlap_CC;
  std::tie(res_proj_solutions, res_overlap_CC) = construct_projected_solutions(solutions, overlap, dims);
  ASSERT_THAT(res_proj_solutions.data(), Pointwise(DoubleEq(), ref_proj_solutions.data()));
  ASSERT_THAT(res_overlap_CC.data(), Pointwise(DoubleEq(), ref_overlap_CC.data()));
}

TEST(itsolv_reset_dspace, select_full_solutions) {
  const size_t nC = 4;
  auto solutions = Matrix<double>({0, 1, 2, 3}, {nC, 1});
  auto norm = std::vector<double>{0.1, 0.4, 0.4, 0.5};
  const double thresh = 0.3;
  const auto ref_full_solutions = std::vector<unsigned int>{1, 2};
  const auto ref_solutions = Matrix<double>({solutions(0, 0), solutions(3, 0)}, {2, 1});
  const auto ref_norm = std::vector<double>{norm[0], norm[3]};
  auto full_solutions = select_full_solutions(solutions, norm, 2, thresh);
  ASSERT_THAT(full_solutions, Pointwise(DoubleEq(), ref_full_solutions));
  ASSERT_THAT(solutions.data(), Pointwise(DoubleEq(), ref_solutions.data()));
  ASSERT_THAT(ref_norm, Pointwise(DoubleEq(), norm));
}