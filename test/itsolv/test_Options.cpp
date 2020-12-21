#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/itsolv/Options.h>

using molpro::linalg::itsolv::Options;

TEST(Options, default_constructor) {
  auto opt = Options{};
  ASSERT_FALSE(opt.n_roots);
  ASSERT_FALSE(opt.convergence_threshold);
}
