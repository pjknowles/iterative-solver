#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/itsolv/LinearEigensystemOptions.h>
#include <molpro/linalg/itsolv/LinearEquationsOptions.h>
#include <molpro/linalg/itsolv/Options.h>

using molpro::linalg::itsolv::LinearEigensystemOptions;
using molpro::linalg::itsolv::LinearEquationsOptions;
using molpro::linalg::itsolv::Options;
using molpro::linalg::itsolv::options_map;

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
  const auto opt_map =
      options_map{{"n_roots", "3"}, {"convergence_threshold", "1.23e-11"}, {"random_test_jhjhaw", "ajwd"}};
  auto opt = Options(opt_map);
  ASSERT_EQ(opt.n_roots.value(), 3);
  ASSERT_DOUBLE_EQ(opt.convergence_threshold.value(), 1.23e-11);
}

TEST(LinearEigensystemOptions, constructor_string) {
  const auto opt_map = options_map{{"n_roots", "3"},
                                   {"convergence_threshold", "1.23e-11"},
                                   {"reset_d", "5"},
                                   {"reset_d_max_q_size", "20"},
                                   {"max_size_qspace", "30"},
                                   {"norm_thresh", "1.3e-4"},
                                   {"svd_thresh", "7.123e-11"},
                                   {"hermiticity", "true"},
                                   {"random_test_12wdgjh", "gi98a"}};
  auto opt = LinearEigensystemOptions(opt_map);
  ASSERT_EQ(opt.n_roots.value(), 3);
  ASSERT_EQ(opt.convergence_threshold.value(), 1.23e-11);
  ASSERT_EQ(opt.reset_D.value(), 5);
  ASSERT_EQ(opt.reset_D_max_Q_size.value(), 20);
  ASSERT_EQ(opt.max_size_qspace.value(), 30);
  ASSERT_DOUBLE_EQ(opt.norm_thresh.value(), 1.3e-4);
  ASSERT_DOUBLE_EQ(opt.svd_thresh.value(), 7.123e-11);
  ASSERT_TRUE(opt.hermiticity.value());
}

TEST(LinearEquationsOptions, constructor_options_map) {
  const auto opt_map = options_map{{"n_roots", "3"},
                                   {"convergence_threshold", "1.23e-11"},
                                   {"reset_d", "5"},
                                   {"reset_d_max_q_size", "20"},
                                   {"max_size_qspace", "30"},
                                   {"norm_thresh", "1.3e-4"},
                                   {"svd_thresh", "7.123e-11"},
                                   {"hermiticity", "true"},
                                   {"augmented_hessian", "0.1"},
                                   {"random_test_12wdgjh", "gi98a"}};
  auto opt = LinearEquationsOptions(opt_map);
  ASSERT_EQ(opt.n_roots.value(), 3);
  ASSERT_EQ(opt.convergence_threshold.value(), 1.23e-11);
  ASSERT_EQ(opt.reset_D.value(), 5);
  ASSERT_EQ(opt.reset_D_max_Q_size.value(), 20);
  ASSERT_EQ(opt.max_size_qspace.value(), 30);
  ASSERT_DOUBLE_EQ(opt.norm_thresh.value(), 1.3e-4);
  ASSERT_DOUBLE_EQ(opt.svd_thresh.value(), 7.123e-11);
  ASSERT_TRUE(opt.hermiticity.value());
  ASSERT_DOUBLE_EQ(opt.augmented_hessian.value(), 0.1);
}
