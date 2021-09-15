#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <molpro/linalg/options.h>

TEST(LinearAlgebra, options) {
  molpro::linalg::set_options(molpro::Options("THINGY", "A=44"));
  const auto options = molpro::linalg::options();
  EXPECT_EQ(options->parameter("A", 99), 44);
  EXPECT_EQ(options->parameter("A", "99"), "44");
}
