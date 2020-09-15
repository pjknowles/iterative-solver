#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/itsolv/subspace/RSpace.h>

using molpro::linalg::itsolv::ArrayHandlers;
using molpro::linalg::itsolv::IterativeSolver;
using molpro::linalg::itsolv::Statistics;
using molpro::linalg::itsolv::subspace::EqnData;
using molpro::linalg::itsolv::subspace::Matrix;
using molpro::linalg::itsolv::subspace::RSpace;
using ::testing::DoubleEq;
using ::testing::Each;
using ::testing::Eq;
using ::testing::Ne;
using ::testing::Pointwise;

/*
 * Things to test:
 *  - update adds new parameters
 *  - update adds new parameters
 */

template <class R, class Q, class P>
struct DummySolver : IterativeSolver<R, Q, P> {
  using typename IterativeSolver<R, Q, P>::value_type;
  using typename IterativeSolver<R, Q, P>::scalar_type;
  size_t add_vector(std::vector<R>& parameters, std::vector<R>& action) override { return 0; };
  size_t add_vector(std::vector<R>& parameters, std::vector<R>& action, std::vector<P>& parametersP) override {
    return 0;
  };
  size_t add_p(std::vector<P>& Pvectors, const value_type* PP, std::vector<R>& parameters, std::vector<R>& action,
               std::vector<P>& parametersP) override {
    return 0;
  };
  void solution(const std::vector<int>& roots, std::vector<R>& parameters, std::vector<R>& residual) override{};
  void solution(const std::vector<int>& roots, std::vector<R>& parameters, std::vector<R>& residual,
                std::vector<P>& parametersP) override{};
  std::vector<size_t> suggest_p(const std::vector<R>& solution, const std::vector<R>& residual, size_t maximumNumber,
                                double threshold) override {
    return {};
  };

  const std::vector<int>& working_set() const override { return ws; };
  const std::vector<scalar_type>& errors() const override { return er; };
  const Statistics& statistics() const override { return st; };
  std::vector<int> ws;
  std::vector<scalar_type> er;
  Statistics st;
};

struct RSpaceF : ::testing::Test {
  using R = std::vector<double>;
  using Q = R;
  using P = std::map<size_t, double>;
  RSpaceF() : handlers(std::make_shared<ArrayHandlers<R, Q, P>>()), rspace(handlers) {}

  std::shared_ptr<ArrayHandlers<R, Q, P>> handlers;
  RSpace<R, Q, P> rspace;
  DummySolver<R, Q, P> solver;
};

TEST_F(RSpaceF, update_null) {
  auto param = std::vector<R>{};
  ASSERT_NO_THROW(rspace.update(param, param, solver));
  ASSERT_EQ(rspace.size(), 0);
  ASSERT_TRUE(rspace.working_set().empty());
  ASSERT_EQ(rspace.data.size(), 2);
  ASSERT_TRUE(rspace.data[EqnData::S].empty());
  ASSERT_TRUE(rspace.data[EqnData::H].empty());
}

TEST_F(RSpaceF, update_single) {
  auto param = std::vector<R>{{1, 2, 3}};
  const auto alpha = 2.0;
  auto action = param;
  for (auto& x : action[0])
    x *= alpha;
  ASSERT_NO_THROW(rspace.update(param, action, solver));
  ASSERT_EQ(rspace.size(), 1);
  ASSERT_THAT(rspace.working_set(), Pointwise(DoubleEq(), std::vector<size_t>{0}));
  const auto& s = rspace.data[EqnData::S];
  const auto& h = rspace.data[EqnData::H];
  ASSERT_EQ(s.size(), 1);
  ASSERT_EQ(h.size(), 1);
  ASSERT_DOUBLE_EQ(s(0, 0), 1.);
  ASSERT_DOUBLE_EQ(h(0, 0), alpha);
}
