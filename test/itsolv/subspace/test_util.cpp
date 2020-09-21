#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/itsolv/subspace/util.h>
#include <numeric>

using molpro::linalg::itsolv::ArrayHandlers;
using molpro::linalg::itsolv::subspace::Matrix;
using molpro::linalg::itsolv::subspace::util::eye_order;
using molpro::linalg::itsolv::subspace::util::overlap;
using molpro::linalg::itsolv::subspace::util::wrap;
using ::testing::DoubleEq;
using ::testing::Each;
using ::testing::Eq;
using ::testing::Ne;
using ::testing::Pointwise;

struct OverlapF : ::testing::Test {
  using R = std::vector<double>;
  using Q = R;
  using P = std::map<size_t, double>;
  OverlapF() : alphas({1, 2, 3}) {
    x.resize(alphas.size());
    for (size_t i = 0; i < alphas.size(); ++i) {
      x[i].assign(nx, alphas[i]);
    }
    ref_overlap.resize({x.size(), x.size()});
    for (size_t i = 0; i < alphas.size(); ++i) {
      for (size_t j = 0; j < alphas.size(); ++j) {
        ref_overlap(i, j) = nx * alphas[i] * alphas[j];
      }
    }
  }
  std::vector<R> x;
  const size_t nx = 5;
  const std::vector<double> alphas;
  Matrix<double> ref_overlap;
  ArrayHandlers<R, Q, P> handler;
};

TEST_F(OverlapF, null_vectors) {
  auto m = overlap<R, R>({}, {}, handler.rr());
  ASSERT_TRUE(m.empty());
}

TEST_F(OverlapF, overlap_one_param) {
  auto m = overlap(wrap(x), handler.rr());
  ASSERT_THAT(m.data(), Pointwise(DoubleEq(), ref_overlap.data()));
}

TEST_F(OverlapF, overlap_two_params) {
  auto m = overlap(wrap(x), wrap(x), handler.rr());
  ASSERT_THAT(m.data(), Pointwise(DoubleEq(), ref_overlap.data()));
}

TEST(EyeOrder, null) {
  auto m = Matrix<double>{};
  auto order = eye_order(m);
  ASSERT_TRUE(order.empty());
}

TEST(EyeOrder, identity) {
  auto&& data = std::vector<double>{1, 0, 0, 0, 1, 0, 0, 0, 1};
  auto m = Matrix<double>{std::move(data), {3, 3}};
  auto order = eye_order(m);
  const auto reference = std::vector<size_t>{0, 1, 2};
  ASSERT_FALSE(order.empty());
  ASSERT_THAT(order, Pointwise(Eq(), reference)) << "matrix is already identity";
}

TEST(EyeOrder, identity_circular_shift) {
  auto&& data = std::vector<double>{0, 1, 0, 0, 0, 1, 1, 0, 0};
  auto m = Matrix<double>{std::move(data), {3, 3}};
  auto order = eye_order(m);
  const auto reference = std::vector<size_t>{2, 0, 1};
  ASSERT_FALSE(order.empty());
  ASSERT_THAT(order, Pointwise(Eq(), reference));
}

TEST(EyeOrder, real_scenario) {
  auto&& data = std::vector<double>{0.1, 0.5, 0.2, 0.2, 0.1, 0.5, 0.5, 0.2, 0.1};
  auto m = Matrix<double>{std::move(data), {3, 3}};
  auto order = eye_order(m);
  const auto reference = std::vector<size_t>{2, 0, 1};
  ASSERT_FALSE(order.empty());
  ASSERT_THAT(order, Pointwise(Eq(), reference));
}