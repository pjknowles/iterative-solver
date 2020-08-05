#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/array/Span.h>

#include <complex>
#include <vector>

using molpro::linalg::array::Span;
using ::testing::Eq;
using ::testing::Pointwise;

template <typename T>
struct SpanF : ::testing::Test {
public:
  SpanF() = default;

  std::vector<T> data = {1, 2, 3, 4, 5};
};

TYPED_TEST_SUITE_P(SpanF);

TYPED_TEST_P(SpanF, default_constructor) {
  auto s = Span<TypeParam>();
  EXPECT_EQ(s.begin(), nullptr);
  EXPECT_EQ(s.size(), 0);
}

TYPED_TEST_P(SpanF, constructor) { auto s = Span<TypeParam>(this->data.data(), this->data.size()); }

TYPED_TEST_P(SpanF, empty) {
  auto& data = this->data;
  auto s = Span<TypeParam>(data.data(), data.size());
  EXPECT_FALSE(s.empty());
  auto t = Span<TypeParam>();
  EXPECT_TRUE(t.empty());
  auto v = Span<TypeParam>(nullptr, 0);
  EXPECT_TRUE(v.empty());
  auto u = Span<TypeParam>(&data[0], 0);
  EXPECT_TRUE(u.empty());
}

TYPED_TEST_P(SpanF, copy_constructor) {
  auto& data = this->data;
  auto s = Span<TypeParam>(data.data(), data.size());
  auto t = Span<TypeParam>(s);
  EXPECT_EQ(s.begin(), t.begin());
  EXPECT_EQ(s.size(), t.size());
}

TYPED_TEST_P(SpanF, move_constructor) {
  auto& data = this->data;
  auto&& s = Span<TypeParam>(data.data(), data.size());
  auto t = Span<TypeParam>(std::forward<Span<TypeParam>>(s));
  EXPECT_EQ(t.begin(), data.data());
  EXPECT_EQ(t.size(), data.size());
  EXPECT_EQ(s.begin(), nullptr);
  EXPECT_EQ(s.size(), 0);
}

TYPED_TEST_P(SpanF, copy_operator) {
  auto& data = this->data;
  auto s = Span<TypeParam>(data.data(), data.size());
  auto t = Span<TypeParam>();
  t = s;
  EXPECT_EQ(s.begin(), t.begin());
  EXPECT_EQ(s.size(), t.size());
}

TYPED_TEST_P(SpanF, move_operator) {
  auto& data = this->data;
  auto&& s = Span<TypeParam>(data.data(), data.size());
  auto t = Span<TypeParam>();
  t = std::forward<Span<TypeParam>>(s);
  EXPECT_EQ(t.begin(), data.data());
  EXPECT_EQ(t.size(), data.size());
  EXPECT_EQ(s.begin(), nullptr);
  EXPECT_EQ(s.size(), 0);
}

TYPED_TEST_P(SpanF, iterate) {
  auto& data = this->data;
  auto ref_data = std::vector<TypeParam>(data.size());
  std::transform(begin(data), end(data), begin(ref_data), [](auto el) { return el + TypeParam(1); });
  auto s = Span<TypeParam>(data.data(), data.size());
  for (auto& el : s)
    el += 1;
  EXPECT_THAT(s, Pointwise(Eq(), ref_data));
  EXPECT_THAT(data, Pointwise(Eq(), ref_data));
  EXPECT_THAT(s, Pointwise(Eq(), data));
}

TYPED_TEST_P(SpanF, begin_end) {
  auto& data = this->data;
  auto ref_data = std::vector<TypeParam>(data.size());
  std::transform(begin(data), end(data), begin(ref_data), [](auto el) { return el + TypeParam(1); });
  auto s = Span<TypeParam>(data.data(), data.size());
  std::transform(begin(s), end(s), begin(s), [](auto el) { return el + TypeParam(1); });
  EXPECT_THAT(s, Pointwise(Eq(), ref_data));
  EXPECT_THAT(data, Pointwise(Eq(), ref_data));
  EXPECT_THAT(s, Pointwise(Eq(), data));
}

REGISTER_TYPED_TEST_SUITE_P(SpanF, default_constructor, constructor, empty, copy_constructor, move_constructor,
                            copy_operator, move_operator, iterate, begin_end);

using ArrayTypes = ::testing::Types<int, double, float, std::complex<double>>;
INSTANTIATE_TYPED_TEST_SUITE_P(Numeric, SpanF, ArrayTypes);
