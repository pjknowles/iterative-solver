#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/itsolv/subspace/RSpace.h>

using molpro::linalg::itsolv::ArrayHandlers;
using molpro::linalg::itsolv::IterativeSolver;
using molpro::linalg::itsolv::Statistics;
using molpro::linalg::itsolv::subspace::EqnData;
using molpro::linalg::itsolv::subspace::Matrix;
using molpro::linalg::itsolv::subspace::RSpace;
using molpro::linalg::itsolv::subspace::rspace::assign_last_parameters_to_new;
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

namespace {
template <class RSpace>
void test_single(RSpace& rspace, double alpha, size_t size, const std::string& message = "") {
  ASSERT_EQ(rspace.size(), size) << message;
  auto ref_working_set = std::vector<size_t>(size);
  std::iota(begin(ref_working_set), end(ref_working_set), size_t{0});
  ASSERT_THAT(rspace.working_set(), Pointwise(DoubleEq(), ref_working_set)) << message;
  const auto& s = rspace.data[EqnData::S];
  const auto& h = rspace.data[EqnData::H];
  ASSERT_EQ(s.size(), size * size) << message;
  ASSERT_EQ(h.size(), size * size) << message;
  for (size_t i = 0; i < size; ++i) {
    ASSERT_DOUBLE_EQ(s(0, 0), 1.) << "i = " << i << " " << message;
    ASSERT_DOUBLE_EQ(h(0, 0), alpha) << "i = " << i << " " << message;
  }
}
} // namespace

TEST_F(RSpaceF, update_single) {
  auto param = std::vector<R>{{1, 2, 3}};
  const auto alpha = 2.0;
  auto action = param;
  for (auto& x : action[0])
    x *= alpha;
  ASSERT_NO_THROW(rspace.update(param, action, solver));
  test_single(rspace, alpha, 1);
}

TEST_F(RSpaceF, update_same_vector_mulitple_times) {
  auto param = std::vector<R>{{1, 2, 3}};
  const auto alpha = 2.0;
  auto action = param;
  for (auto& x : action[0])
    x *= alpha;
  for (size_t i = 0; i < 4; ++i) {
    ASSERT_NO_THROW(rspace.update(param, action, solver));
    test_single(rspace, alpha, 1);
    rspace.update_working_set({0});
  }
}

TEST_F(RSpaceF, assign_new_parameters_to_last__stable_ordering) {
  auto last_param = std::vector<R>{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  auto assign_param_to_last = assign_last_parameters_to_new(last_param, last_param, handlers->qr());
  ASSERT_THAT(assign_param_to_last, Pointwise(Eq(), std::vector<size_t>{0, 1, 2}))
      << "same params should be in ascending order";
  auto param = last_param;
  std::swap(param[0], param[2]);
  assign_param_to_last = assign_last_parameters_to_new(param, last_param, handlers->qr());
  ASSERT_THAT(assign_param_to_last, Pointwise(Eq(), std::vector<size_t>{2, 1, 0})) << " reversed parameters";
  param = last_param;
  std::swap(param[0], param[1]);
  std::swap(param[1], param[2]);
  assign_param_to_last = assign_last_parameters_to_new(param, last_param, handlers->qr());
  ASSERT_THAT(assign_param_to_last, Pointwise(Eq(), std::vector<size_t>{1, 2, 0})) << " cyclic permutation";
}

TEST_F(RSpaceF, update_working_set) {
  auto param = std::vector<R>{{1, 2, 3}};
  const auto alpha = 2.0;
  auto action = param;
  for (auto& x : action[0])
    x *= alpha;
  ASSERT_NO_THROW(rspace.update(param, action, solver));
  test_single(rspace, alpha, 1);
  rspace.update_working_set({0});
  test_single(rspace, alpha, 1);
  ASSERT_EQ(rspace.last_params().size(), 1);
  ASSERT_EQ(rspace.last_actions().size(), 1);
  rspace.update_working_set({});
  ASSERT_EQ(rspace.last_params().size(), 0);
  ASSERT_EQ(rspace.last_actions().size(), 0);
  ASSERT_EQ(rspace.working_set().size(), 0);
}

TEST_F(RSpaceF, update_check_ordering) {
  auto param = std::vector<R>{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  const auto size = param.size();
  ASSERT_NO_THROW(rspace.update(param, param, solver));
  test_single(rspace, 1, size);
  rspace.update_working_set({0, 1, 2});
  ASSERT_NO_THROW(rspace.update(param, param, solver));
  test_single(rspace, 1, size, "update with the same parameters should leave order unchanged");
  rspace.update_working_set({0, 1, 2});
  auto param_reverse = param;
  std::swap(param_reverse[0], param_reverse[2]);
  ASSERT_NO_THROW(rspace.update(param_reverse, param_reverse, solver));
  ASSERT_EQ(rspace.size(), size);
  ASSERT_THAT(rspace.working_set(), Pointwise(DoubleEq(), std::vector<size_t>{0, 1, 2}))
      << "changing order of input parameters does not change the order of the working set";
}
