#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "DummySolver.h"
#include <molpro/linalg/itsolv/subspace/QSpace.h>

using molpro::linalg::itsolv::ArrayHandlers;
using molpro::linalg::itsolv::subspace::EqnData;
using molpro::linalg::itsolv::subspace::Matrix;
using molpro::linalg::itsolv::subspace::null_data;
using molpro::linalg::itsolv::subspace::QSpace;
using molpro::linalg::itsolv::subspace::RSpace;
using molpro::linalg::itsolv::subspace::SubspaceData;
using molpro::linalg::itsolv::subspace::qspace::merge_subspace_qq;
using molpro::linalg::itsolv::subspace::qspace::merge_subspace_qr;
using molpro::linalg::itsolv::subspace::qspace::merge_subspace_rq;
using molpro::linalg::itsolv::subspace::qspace::QParam;
using molpro::linalg::itsolv::subspace::qspace::update;
using molpro::linalg::itsolv::subspace::qspace::update_qq_subspace;
using molpro::linalg::itsolv::subspace::qspace::update_qr_subspace;
using molpro::linalg::itsolv::subspace::util::wrap;
using ::testing::DoubleEq;
using ::testing::Each;
using ::testing::Eq;
using ::testing::Ne;
using ::testing::Pointwise;

namespace {
using R = std::vector<double>;
using Q = R;
using P = std::map<size_t, double>;
} // namespace

// Simple test case
struct QSpaceUpdateF : ::testing::Test {
  auto call_update() {
    if (qparam.empty() && !params.empty()) {
      qparam.assign(params[0].size(), 1);
      qaction = qparam;
    }
    if (working_set.empty()) {
      working_set.resize(params.size());
      std::iota(begin(working_set), end(working_set), size_t{0});
    }
    return update(qparam, qaction, wrap(params), wrap(actions), last_params, last_actions, working_set, handlers);
  }

  Q qparam{};
  Q qaction{};
  std::vector<R> params;
  std::vector<R> actions;
  std::vector<R> last_params;
  std::vector<R> last_actions;
  std::vector<size_t> working_set;
  DummySolver<R, Q, P> solver;
  ArrayHandlers<R, Q, P> handlers;
};

TEST_F(QSpaceUpdateF, update_null) {
  auto result = call_update();
  auto &result_params = result.first;
  auto &used_working_set = result.second;
  ASSERT_TRUE(result_params.empty());
  ASSERT_TRUE(used_working_set.empty());
}

TEST_F(QSpaceUpdateF, update_single_rspace_zero) {
  auto zero = std::vector<R>{{0, 0, 0, 0}};
  params = zero;
  actions = zero;
  last_params = zero;
  last_actions = zero;
  auto result = call_update();
  auto &result_params = result.first;
  auto &used_working_set = result.second;
  ASSERT_TRUE(result_params.empty());
  ASSERT_TRUE(used_working_set.empty());
}

TEST_F(QSpaceUpdateF, update_single_rspace_equal_new_and_last) {
  auto one = std::vector<R>{{1, 0, 0, 0}};
  params = one;
  actions = one;
  last_params = one;
  last_actions = one;
  auto result = call_update();
  auto &result_params = result.first;
  auto &used_working_set = result.second;
  ASSERT_TRUE(result_params.empty()) << "difference vector should be zero and must not be added";
  ASSERT_TRUE(used_working_set.empty());
}

TEST_F(QSpaceUpdateF, update_single_rspace_zero_last) {
  auto zero = std::vector<R>{{0, 0, 0, 0}};
  auto one = std::vector<R>{{1, 0, 0, 0}};
  params = one;
  actions = one;
  last_params = zero;
  last_actions = zero;
  auto result = call_update();
  auto &result_params = result.first;
  auto &used_working_set = result.second;
  ASSERT_EQ(result_params.size(), 1)
      << "new params are valid and last params are zero; should lead to a single difference vector";
  ASSERT_EQ(used_working_set.size(), 1);
}

TEST_F(QSpaceUpdateF, update_multiple_rspace_zero_last) {
  auto zero = std::vector<R>{{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
  auto eye = std::vector<R>{{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
  params = eye;
  actions = eye;
  last_params = zero;
  last_actions = zero;
  auto result = call_update();
  auto &result_params = result.first;
  auto &used_working_set = result.second;
  ASSERT_EQ(result_params.size(), params.size())
      << "new params are valid and last params are zero; each param should give a valid difference vector";
  ASSERT_EQ(used_working_set.size(), params.size());
}

TEST_F(QSpaceUpdateF, update_difference_reproduces_new_params) {
  params = {{1, 2, 3}};
  actions = params;
  last_params = {{3, 2, 1}};
  last_actions = last_params;
  auto result = call_update();
  auto &result_params = result.first;
  auto &used_working_set = result.second;
  ASSERT_EQ(result_params.size(), params.size());
  ASSERT_EQ(used_working_set.size(), params.size());
  auto &qdiff_params = *result_params.front().param;
  auto norm = result_params.front().normalisation_constant;
  handlers.rr().scal(1. / norm, qdiff_params);
  handlers.rr().axpy(1., last_params.front(), qdiff_params);
  ASSERT_THAT(qdiff_params, Pointwise(DoubleEq(), params.front()))
      << "adding difference vector to last should give back new params";
}

TEST(QSpaceUpdateSubspace, qq) {
  auto old_params = std::vector<Q>{{10}, {20}};
  auto new_params = std::vector<Q>{{1}, {2}, {3}};
  SubspaceData qq = null_data<EqnData::H, EqnData::S>();
  //@formatter:off
  auto reference = std::vector<double>{
     0,  0, 10, 20, 30,
     0,  0, 20, 40, 60,
    10, 20,  1,  2,  3,
    20, 40,  2,  4,  6,
    30, 60,  3,  6,  9
  };
  //@formatter:on
  auto handler = ArrayHandlers<R, Q, P>{};
  update_qq_subspace(wrap(old_params), wrap(old_params), wrap(new_params), wrap(new_params), qq, handler.qq());
  ASSERT_THAT(qq[EqnData::S].data(), Pointwise(DoubleEq(), reference));
  ASSERT_THAT(qq[EqnData::H].data(), Pointwise(DoubleEq(), reference));
}

TEST(QSpaceUpdateSubspace, qr) {
  auto rparams = std::vector<R>{{10}, {20}};
  auto qparams = std::vector<Q>{{1}, {2}, {3}};
  SubspaceData qr = null_data<EqnData::H, EqnData::S>();
  SubspaceData rq = null_data<EqnData::H, EqnData::S>();
  //@formatter:off
  auto reference_qr = std::vector<double>{
      10, 20,
      20, 40,
      30, 60
  };
  auto reference_rq = std::vector<double>{
      10, 20, 30,
      20, 40, 60
  };
  //@formatter:on
  auto handler = ArrayHandlers<R, Q, P>{};
  update_qr_subspace(wrap(qparams), wrap(qparams), wrap(rparams), wrap(rparams), qr, rq, handler.qr(), handler.rq());
  ASSERT_THAT(qr[EqnData::S].data(), Pointwise(DoubleEq(), reference_qr));
  ASSERT_THAT(qr[EqnData::H].data(), Pointwise(DoubleEq(), reference_qr));
  ASSERT_THAT(rq[EqnData::S].data(), Pointwise(DoubleEq(), reference_rq));
  ASSERT_THAT(rq[EqnData::H].data(), Pointwise(DoubleEq(), reference_rq));
}

//this is an overkill, but might be useful for other tests
struct QSpaceEyeF : ::testing::Test {
  QSpaceEyeF() : qspace(handlers) {}

  void SetUp() override {
    auto zero = R{0, 0, 0};
    auto eye = std::vector<R>{{1, 0, 0, 0, 0}, {0, 1, 0, 0, 0}, {0, 0, 1, 0, 0}, {0, 0, 0, 1, 0}, {0, 0, 0, 0, 1}};
    auto init_params = std::vector<R>{eye[0], eye[0]};
    RSpace<R, Q, P> rspace{handlers};
    rspace.update(init_params, init_params, solver);
    rspace.update_working_set({0, 1});
    auto rparams = std::list<std::vector<R>>{
        {eye[1], eye[1]},
        {eye[2], eye[2]},
        {eye[3], eye[3]},
    };
    for (auto &p : rparams) {
      rspace.update(p, p, solver);
      qspace.update(rspace, solver);
      rspace.update_working_set({0, 1});
      ASSERT_EQ(qspace.used_working_set().size(), 2) << "building up q space two roots at a time";
    }
    auto rparams_root_1 = std::vector<R>{eye[4]};
    rspace.update_working_set({1});
    rspace.update(rparams_root_1, rparams_root_1, solver);
    qspace.update(rspace, solver);
    rspace.update_working_set({0});
    ASSERT_EQ(qspace.used_working_set().size(), 1) << "only root 1 at the end";
  }

  std::shared_ptr<ArrayHandlers<R, Q, P>> handlers = std::make_shared<ArrayHandlers<R, Q, P>>();
  QSpace<R, Q, P> qspace;
  std::vector<size_t> working_set = {0, 1};
  DummySolver<R, Q, P> solver;
};

TEST_F(QSpaceEyeF, modification_candidates_null) {
  auto qs = QSpace<R, Q, P>{handlers};
  auto candidates = qs.modification_candidates(0);
  ASSERT_TRUE(candidates.empty());
}

TEST_F(QSpaceEyeF, modification_candidates) {
  auto candidates = qspace.modification_candidates(0);
  ASSERT_FALSE(candidates.empty());
  auto reference = std::vector<size_t>{0, 2};
  ASSERT_THAT(candidates, Pointwise(Eq(), reference));
  reference = std::vector<size_t>{1, 3, 5};
  candidates = qspace.modification_candidates(1);
  ASSERT_THAT(candidates, Pointwise(Eq(), reference));
}

struct QSpaceMergeF : ::testing::Test {
  QSpaceMergeF() {
    auto dummy = QParam<Q>{};
    for (auto d : {EqnData::H, EqnData::S}) {
      qq[d].resize({nq, nq});
      qr[d].resize({nq, nr});
      rq[d].resize({nr, nq});
      qq[d].fill(alpha);
      qr[d].fill(alpha);
      rq[d].fill(alpha);
    }
  }

  SubspaceData qq =null_data<EqnData::H,EqnData::S>();
  SubspaceData qr =null_data<EqnData::H,EqnData::S>();
  SubspaceData rq =null_data<EqnData::H,EqnData::S>();
  const size_t nq = 3;
  const size_t nr = 2;
  const double alpha = 1.;
  const double a = 2.;
  const double b = 3.;
};

TEST_F(QSpaceMergeF, merge_qq) {
  const size_t i = 0, j = 1;
  merge_subspace_qq(i, j, a, b, qq);
  auto dim_ref = Matrix<double>::coord_type{nq - 1, nq - 1};
  auto ref_mat = Matrix<double>(dim_ref);
  ref_mat.fill(alpha);
  const auto result = a * alpha + b * alpha;
  const auto diag = result * result;
  ref_mat.slice({i, 0}, {i + 1, ref_mat.cols()}).fill(result);
  ref_mat.slice({0, i}, {ref_mat.rows(), i + 1}).fill(result);
  ref_mat(i,i) = diag;
  for (auto d : {EqnData::H, EqnData::S}) {
    ASSERT_EQ(qq[d].dimensions(), dim_ref);
    ASSERT_THAT(qq[d].data(), Pointwise(DoubleEq(), ref_mat.data()));
  }
}

TEST_F(QSpaceMergeF, merge_qr) {
  const size_t i = 0, j = 1;
  merge_subspace_qr(i, j, a, b, qr);
  auto dim_ref = Matrix<double>::coord_type{nq - 1, nr};
  auto ref_mat = Matrix<double>(dim_ref);
  const auto result = a * alpha + b * alpha;
  ref_mat.fill(alpha);
  ref_mat.slice({i, 0}, {i + 1, ref_mat.cols()}).fill(result);
  for (auto d : {EqnData::H, EqnData::S}) {
    ASSERT_EQ(qr[d].dimensions(), dim_ref);
    ASSERT_THAT(qr[d].data(), Pointwise(DoubleEq(), ref_mat.data()));
  }
}

TEST_F(QSpaceMergeF, merge_rq) {
  const size_t i = 0, j = 1;
  merge_subspace_rq(i, j, a, b, rq);
  auto dim_ref = Matrix<double>::coord_type{nr, nq - 1};
  auto ref_mat = Matrix<double>(dim_ref);
  const auto result = a * alpha + b * alpha;
  ref_mat.fill(alpha);
  ref_mat.slice({0, i}, {ref_mat.rows(), i + 1}).fill(result);
  for (auto d : {EqnData::H, EqnData::S}) {
    ASSERT_EQ(rq[d].dimensions(), dim_ref);
    ASSERT_THAT(rq[d].data(), Pointwise(DoubleEq(), ref_mat.data()));
  }
}
