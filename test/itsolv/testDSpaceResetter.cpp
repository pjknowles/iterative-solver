#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "subspace/DummyXSpace.h"
#include <molpro/linalg/itsolv/DSpaceResetter.h>

using molpro::linalg::itsolv::cwrap;
using molpro::linalg::itsolv::wrap;
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
