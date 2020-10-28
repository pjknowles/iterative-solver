#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/itsolv/reset_dspace.h>

using molpro::linalg::itsolv::detail::construct_projected_solutions;
using molpro::linalg::itsolv::detail::DoReset;
using molpro::linalg::itsolv::subspace::Matrix;
using molpro::linalg::itsolv::subspace::xspace::Dimensions;
using ::testing::DoubleEq;
using ::testing::Pointwise;

struct DoResetF : public ::testing::Test, DoReset {};

TEST(itsolv_reset_dspace, DoReset_constructor) {
  auto d = DoReset();
  ASSERT_FALSE(d.value());
  const size_t n = 7;
  auto d2 = DoReset(n);
  ASSERT_FALSE(d2.value());
}

TEST(itsolv_reset_dspace, DoReset_set_get_nreset) {
  const size_t n = 7, m = 3;
  auto d = DoReset(n);
  ASSERT_EQ(d.get_nreset(), n);
  d.set_nreset(m);
  ASSERT_EQ(d.get_nreset(), m);
}

TEST(itsolv_reset_dspace, DoReset_update) {
  const size_t n = 2;
  const unsigned int init_max_q = 10, nD = 3;
  unsigned int max_q = init_max_q;
  auto d = DoReset(n);
  d.update(0, max_q, nD);
  ASSERT_FALSE(d.value());
  ASSERT_EQ(max_q, init_max_q);
  d.update(1, max_q, nD);
  ASSERT_TRUE(d.value());
  ASSERT_EQ(max_q, init_max_q + nD);
  d.update(2, max_q, nD);
  ASSERT_TRUE(d.value());
  ASSERT_EQ(max_q, init_max_q + nD);
  d.update(3, max_q, 0);
  ASSERT_FALSE(d.value());
  ASSERT_EQ(max_q, init_max_q);
}

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
