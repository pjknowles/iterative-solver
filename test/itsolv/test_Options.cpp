#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/itsolv/Options.h>

using molpro::linalg::itsolv::Options;

TEST(Options, default_constructor) {
  auto opt = Options{};
  ASSERT_FALSE(opt.n_roots);
  ASSERT_FALSE(opt.convergence_threshold);
}

TEST(Options, copy) {
  const int n = 3;
  const double conv_thresh = 1.5e-15;
  auto opt = Options{};
  opt.n_roots = n;
  opt.convergence_threshold = conv_thresh;
  auto opt_copy = Options();
  opt_copy.copy(opt);
  ASSERT_EQ(opt_copy.n_roots, opt.n_roots);
  ASSERT_EQ(opt_copy.convergence_threshold, opt.convergence_threshold);
}

TEST(Options, constructor_string) {
  const auto options_map = std::map<std::string, std::string>{{"n_roots", "3"}, {"convergence_threshold", "1.23e-11"}};
  auto opt = Options(options_map);
  ASSERT_EQ(opt.n_roots.value(), 3);
  ASSERT_EQ(opt.convergence_threshold.value(), 1.23e-11);
}
