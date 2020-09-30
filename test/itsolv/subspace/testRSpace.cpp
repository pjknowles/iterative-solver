#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "DummySolver.h"
#include <molpro/linalg/itsolv/subspace/RSpace.h>

using molpro::linalg::itsolv::ArrayHandlers;
using molpro::linalg::itsolv::Logger;
using molpro::linalg::itsolv::subspace::EqnData;
using molpro::linalg::itsolv::subspace::Matrix;
using molpro::linalg::itsolv::subspace::RSpace;
using ::testing::DoubleEq;
using ::testing::Each;
using ::testing::Eq;
using ::testing::Ne;
using ::testing::Pointwise;

struct RSpaceF : ::testing::Test {
  using R = std::vector<double>;
  using Q = R;
  using P = std::map<size_t, double>;
  RSpaceF() : handlers(std::make_shared<ArrayHandlers<R, Q, P>>()), rspace(handlers, logger) {}

  std::shared_ptr<ArrayHandlers<R, Q, P>> handlers;
  std::shared_ptr<Logger> logger = std::make_shared<Logger>();
  RSpace<R, Q, P> rspace;
  DummySolver<R, Q, P> solver;
};

TEST_F(RSpaceF, update_null) {
  auto param = std::vector<R>{};
  ASSERT_NO_THROW(rspace.update(param, param, solver));
  ASSERT_EQ(rspace.size(), 0);
  ASSERT_EQ(rspace.data.size(), 2);
  ASSERT_TRUE(rspace.data[EqnData::S].empty());
  ASSERT_TRUE(rspace.data[EqnData::H].empty());
}

namespace {
template <class RSpace>
void test_single(RSpace& rspace, double alpha, size_t size, size_t working_set_size, const std::string& message = "") {
  ASSERT_EQ(rspace.size(), size) << message;
  auto ref_working_set = std::vector<size_t>(working_set_size);
  std::iota(begin(ref_working_set), end(ref_working_set), size_t{0});
  const auto& s = rspace.data[EqnData::S];
  const auto& h = rspace.data[EqnData::H];
  ASSERT_EQ(s.size(), size * size) << message;
  ASSERT_EQ(h.size(), size * size) << message;
  for (size_t i = 0; i < size; ++i) {
    ASSERT_DOUBLE_EQ(s(i, i), 1.) << "i = " << i << " " << message;
    ASSERT_DOUBLE_EQ(h(i, i), alpha) << "i = " << i << " " << message;
  }
}
} // namespace

TEST_F(RSpaceF, update_single) {
  const double norm = std::sqrt(1. + 4. + 9.);
  auto param = std::vector<R>{{1. / norm, 2. / norm, 3. / norm}};
  const auto alpha = 2.0;
  auto action = param;
  for (auto& x : action[0])
    x *= alpha;
  ASSERT_NO_THROW(rspace.update(param, action, solver));
  test_single(rspace, alpha, 0, 1);
}

TEST_F(RSpaceF, update_same_vector_mulitple_times) {
  const double norm = std::sqrt(1. + 4. + 9.);
  auto param = std::vector<R>{{1. / norm, 2. / norm, 3. / norm}};
  const auto alpha = 2.0;
  auto action = param;
  for (auto& x : action[0])
    x *= alpha;
  for (size_t i = 0; i < 4; ++i) {
    ASSERT_NO_THROW(rspace.update(param, action, solver));
    if (i == 0)
      test_single(rspace, alpha, 0, 1);
    else
      test_single(rspace, alpha, 1, 1);
  }
}
