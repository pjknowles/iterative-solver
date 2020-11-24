#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "subspace/DummyXSpace.h"
#include <molpro/linalg/itsolv/DSpaceResetter.h>

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
