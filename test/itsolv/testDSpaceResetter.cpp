#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "subspace/DummyXSpace.h"
#include <molpro/linalg/itsolv/DSpaceResetter.h>

#include <cmath>

using molpro::linalg::itsolv::cwrap;
using molpro::linalg::itsolv::wrap;
using molpro::linalg::itsolv::detail::DSpaceResetter;
using molpro::linalg::itsolv::detail::max_overlap_with_R;
using molpro::linalg::itsolv::detail::resize_qspace;
using molpro::linalg::itsolv::subspace::Matrix;

using R = std::vector<double>;
using Q = std::vector<double>;
using P = std::map<size_t, double>;

struct ResizeQspace : DummyXSpace<R, Q, P> {
  ResizeQspace() { dims = molpro::linalg::itsolv::subspace::Dimensions{0, static_cast<size_t>(nQ), 0}; }
  void eraseq(size_t i) override {
    --nQ;
    erased_Q.emplace_back(i);
    dims = molpro::linalg::itsolv::subspace::Dimensions(0, nQ, 0);
  }

  int nQ = 3;
  std::vector<size_t> erased_Q;
};

TEST(DSpaceResetter, resize_qspace) {
  auto solutions = Matrix<double>{{0.0, 0.5, -0.3, //
                                   0.0, 0.1, -0.4, //
                                   0.1, 0.2, 0.1},
                                  {3, 3}};
  const size_t max_Qsize_after_reset = 0;
  molpro::linalg::itsolv::Logger logger{};
  auto xspace = ResizeQspace{};
  resize_qspace(xspace, solutions, max_Qsize_after_reset, logger);
  const auto reference_erased_Q = std::vector<size_t>{2, 1, 0};
  ASSERT_EQ(xspace.nQ, max_Qsize_after_reset);
  ASSERT_THAT(xspace.erased_Q, ::testing::Pointwise(::testing::Eq(), reference_erased_Q));
}

TEST(DSpaceResetter, do_reset_periodic) {
  auto resetter = DSpaceResetter<Q>{3, 10};
  const auto dims = molpro::linalg::itsolv::subspace::Dimensions{0, 4, 2};
  EXPECT_FALSE(resetter.do_reset(0, dims));
  EXPECT_FALSE(resetter.do_reset(1, dims));
  EXPECT_TRUE(resetter.do_reset(2, dims));
  EXPECT_FALSE(resetter.do_reset(2, molpro::linalg::itsolv::subspace::Dimensions{0, 4, 0}));
}

TEST(DSpaceResetter, do_reset_never_if_nreset_zero) {
  auto resetter = DSpaceResetter<Q>{0, 10};
  const auto dims = molpro::linalg::itsolv::subspace::Dimensions{0, 4, 2};
  for (size_t iter = 0; iter < 10; ++iter)
    EXPECT_FALSE(resetter.do_reset(iter, dims));
}

// Orthonormal Q and D vectors, whose Q space is reduced after the reset
struct ResetXSpace : DummyXSpace<R, Q, P> {
  ResetXSpace() {
    for (size_t i = 0; i < 4; ++i) {
      auto e = R(4, 0.);
      e[i] = 1;
      (i < 3 ? q : d).push_back(e);
    }
    update_dims();
    data[molpro::linalg::itsolv::subspace::EqnData::S] = Matrix<double>{{1, 0, 0, 0, //
                                                                         0, 1, 0, 0, //
                                                                         0, 0, 1, 0, //
                                                                         0, 0, 0, 1},
                                                                        {4, 4}};
  }
  void update_dims() { dims = molpro::linalg::itsolv::subspace::Dimensions(0, q.size(), d.size()); }
  void eraseq(size_t i) override {
    erased_Q.push_back(i);
    q.erase(q.begin() + i);
    update_dims();
  }
  void update_dspace(VecRef<Q>& params, VecRef<Q>& actions) override {
    d.clear();
    for (const auto& p : params)
      d.push_back(p.get());
    update_dims();
  }
  CVecRef<Q> cparamsq() const override { return cwrap(q); }
  CVecRef<Q> cparamsd() const override { return cwrap(d); }

  std::vector<Q> q, d;
  std::vector<size_t> erased_Q;
};

// The Q vector with the largest overlap with the solution is replaced by it, and the remaining Q space has to be
// reduced by removing the vector with the smallest contribution to the solution, which must be identified in the
// subspace as it is after the erasure and with the D space cleared.
TEST(DSpaceResetter, run_resizes_qspace_after_erasing) {
  auto xspace = ResetXSpace{};
  const auto solutions = Matrix<double>{{0.9, 0.05, 0.3, 0.3}, {1, 4}};
  auto handlers = molpro::linalg::itsolv::ArrayHandlers<R, Q, P>{};
  molpro::linalg::itsolv::Logger logger{};
  auto rparams = std::vector<R>(1, R(4, 0.));
  auto resetter = DSpaceResetter<Q>{1, 2};
  ASSERT_TRUE(resetter.do_reset(0, xspace.dimensions()));
  const auto working_set = resetter.run(wrap(rparams), xspace, solutions, 1.0e-14, 1.0e-12, handlers, logger);
  EXPECT_THAT(working_set, ::testing::ElementsAre(0));
  EXPECT_THAT(xspace.erased_Q, ::testing::ElementsAre(0, 0));
  ASSERT_EQ(xspace.q.size(), 1);
  EXPECT_THAT(xspace.q.front(), ::testing::ElementsAre(0, 0, 1, 0));
  EXPECT_TRUE(xspace.d.empty());
  const double norm = std::sqrt(0.9 * 0.9 + 0.05 * 0.05 + 0.3 * 0.3 + 0.3 * 0.3);
  EXPECT_THAT(rparams.front(),
              ::testing::Pointwise(::testing::DoubleNear(1e-14),
                                   std::vector<double>{0.9 / norm, 0.05 / norm, 0.3 / norm, 0.3 / norm}));
}

struct DSpaceResetter_MaxOverlapWithR_TrivialF : ::testing::Test {
  std::shared_ptr<molpro::linalg::array::ArrayHandler<R, Q>> handler =
      molpro::linalg::array::create_default_handler<R, Q>();
  molpro::linalg::itsolv::Logger logger{};
  std::vector<R> rparams{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
};

TEST_F(DSpaceResetter_MaxOverlapWithR_TrivialF, qparams_3) {
  auto qparams = std::vector<R>{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  auto max_overlap_indices_ref = std::vector<int>{2, 1, 0};
  auto max_overlap_indices = max_overlap_with_R(cwrap(rparams), cwrap(qparams), *handler, logger);
  ASSERT_THAT(max_overlap_indices, ::testing::Pointwise(::testing::Eq(), max_overlap_indices_ref));
}

TEST_F(DSpaceResetter_MaxOverlapWithR_TrivialF, qparams_2) {
  auto qparams = std::vector<R>{{1, 0, 0}, {0, 0, 1}};
  auto max_overlap_indices_ref = std::vector<int>{1, 0};
  auto max_overlap_indices = max_overlap_with_R(cwrap(rparams), cwrap(qparams), *handler, logger);
  ASSERT_THAT(max_overlap_indices, ::testing::Pointwise(::testing::Eq(), max_overlap_indices_ref));
}

TEST_F(DSpaceResetter_MaxOverlapWithR_TrivialF, qparams_1) {
  auto qparams = std::vector<R>{{0, 1, 0}};
  auto max_overlap_indices_ref = std::vector<int>{0};
  auto max_overlap_indices = max_overlap_with_R(cwrap(rparams), cwrap(qparams), *handler, logger);
  ASSERT_THAT(max_overlap_indices, ::testing::Pointwise(::testing::Eq(), max_overlap_indices_ref));
}

TEST_F(DSpaceResetter_MaxOverlapWithR_TrivialF, qparams_0) {
  auto qparams = std::vector<R>{};
  auto max_overlap_indices_ref = std::vector<int>{};
  auto max_overlap_indices = max_overlap_with_R(cwrap(rparams), cwrap(qparams), *handler, logger);
  ASSERT_THAT(max_overlap_indices, ::testing::Pointwise(::testing::Eq(), max_overlap_indices_ref));
}
