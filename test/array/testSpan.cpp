#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/array/Span.h>

#include <vector>

using molpro::linalg::array::Span;
using ::testing::Eq;
using ::testing::Pointwise;

struct SpanF : ::testing::Test {
public:
  SpanF() = default;

  using T = int;
  std::vector<T> data = {1, 2, 3, 4, 5};
};

TEST_F(SpanF, default_constructor) {
  auto s = Span<T>();
  EXPECT_EQ(s.begin(), nullptr);
  EXPECT_EQ(s.size(), 0);
}

TEST_F(SpanF, constructor) { auto s = Span<T>(data.data(), data.size()); }

TEST_F(SpanF, copy_constructor) {
  auto s = Span<T>(data.data(), data.size());
  auto t = Span<T>(s);
  EXPECT_EQ(s.begin(), t.begin());
  EXPECT_EQ(s.size(), t.size());
}

TEST_F(SpanF, move_constructor) {
  auto&& s = Span<T>(data.data(), data.size());
  auto t = Span<T>(std::forward<Span<T>>(s));
  EXPECT_EQ(t.begin(), data.data());
  EXPECT_EQ(t.size(), data.size());
  EXPECT_EQ(s.begin(), nullptr);
  EXPECT_EQ(s.size(), 0);
}

TEST_F(SpanF, copy_operator) {
  auto s = Span<T>(data.data(), data.size());
  auto t = Span<T>();
  t = s;
  EXPECT_EQ(s.begin(), t.begin());
  EXPECT_EQ(s.size(), t.size());
}

TEST_F(SpanF, move_operator) {
  auto&& s = Span<T>(data.data(), data.size());

  auto t = Span<T>();
  t = std::forward<Span<T>>(s);
  EXPECT_EQ(t.begin(), data.data());
  EXPECT_EQ(t.size(), data.size());
  EXPECT_EQ(s.begin(), nullptr);
  EXPECT_EQ(s.size(), 0);
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
