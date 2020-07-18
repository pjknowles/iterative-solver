#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/array/Span.h>

#include <vector>

using molpro::linalg::array::Span;
using ::testing::Eq;
using ::testing::Pointwise;

TEST(Span, constructor) {
  using T = int;
  auto data = std::vector<T>{1, 2, 3, 4, 5};
  auto s = Span<T>(data.data(), data.size());
}

TEST(Span, iterate) {
  using T = int;
  auto data = std::vector<T>{1, 2, 3, 4, 5};
  auto ref_data = std::vector<T>{2, 3, 4, 5, 6};
  auto s = Span<T>(data.data(), data.size());
  for (auto& el : s)
    el += 1;
  EXPECT_THAT(s, Pointwise(Eq(), ref_data));
  EXPECT_THAT(data, Pointwise(Eq(), ref_data));
  EXPECT_THAT(s, Pointwise(Eq(), data));
}
