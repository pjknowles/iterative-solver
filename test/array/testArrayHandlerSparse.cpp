#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/array/ArrayHandlerSparse.h>

using molpro::linalg::array::ArrayHandlerSparse;

using ::testing::ContainerEq;
using ::testing::DoubleEq;
using ::testing::Pointwise;

template <class T>
struct TestArrayHandlerSparse : ::testing::Test {};

TYPED_TEST_SUITE_P(TestArrayHandlerSparse);

TYPED_TEST_P(TestArrayHandlerSparse, constructor) {
  auto handler = ArrayHandlerSparse<std::map<size_t, typename TypeParam::first_type>,
                                    std::map<size_t, typename TypeParam::second_type>>{};
}

TYPED_TEST_P(TestArrayHandlerSparse, dot) {
  auto x = std::map<size_t, typename TypeParam::first_type>{{1, 1.0}, {3, 2.0}, {5, 3.0}, {11, 4.0}};
  auto y = std::map<size_t, typename TypeParam::second_type>{{1, 1.0}, {3, 2.0}, {6, 3.0}, {11, 4.0}, {20, 5.0}};
  auto handler = ArrayHandlerSparse<decltype(x), decltype(y)>{};
  auto ref = typename decltype(handler)::value_type(1. + 2. * 2. + 4. * 4.);
  auto result = handler.dot(x, y);
  ASSERT_DOUBLE_EQ(result, ref);
}

TYPED_TEST_P(TestArrayHandlerSparse, select_max_dot) {
  auto x = std::map<size_t, typename TypeParam::first_type>{{1, -2}, {3, 1}, {4, 3}, {6, -4}};
  auto y = std::map<size_t, typename TypeParam::second_type>{{0, 1}, {1, 1}, {2, 1}, {4, 1}, {6, 1}};
  auto handler = ArrayHandlerSparse<decltype(x), decltype(y)>{};
  auto ref_result = std::map<size_t, typename decltype(handler)::value_type_abs>{{6, 4}, {4, 3}, {1, 2}};
  auto select = handler.select_max_dot(ref_result.size(), x, y);
  ASSERT_THAT(select, ContainerEq(ref_result));
}

TYPED_TEST_P(TestArrayHandlerSparse, axpy) {
  const int alpha = 3;
  auto x = std::map<size_t, typename TypeParam::first_type>{{1, 1.0}, {3, 2.0}, {5, 3.0}, {11, 4.0}};
  auto y = std::map<size_t, typename TypeParam::second_type>{{1, 1.0}, {3, 2.0}, {6, 3.0}, {11, 4.0}, {20, 5.0}};
  auto y_ref = std::map<size_t, typename TypeParam::second_type>{
      {1, y[1] + x[1] * alpha}, {3, y[3] + x[3] * alpha}, {6, y[6]}, {11, y[11] + x[11] * alpha}, {20, y[20]}};
  auto handler = ArrayHandlerSparse<decltype(y), decltype(x)>{};
  handler.axpy(alpha, x, y);
  ASSERT_EQ(y.size(), y_ref.size());
  auto iy = y.begin();
  auto iy_ref = y_ref.begin();
  for (auto i = 0; i < y.size(); ++i, ++iy, ++iy_ref) {
    EXPECT_DOUBLE_EQ(iy->second, iy_ref->second) << "i = " << std::to_string(i);
    EXPECT_EQ(iy->first, iy_ref->first) << "i = " << std::to_string(i);
  }
}

REGISTER_TYPED_TEST_SUITE_P(TestArrayHandlerSparse, constructor, dot, axpy, select_max_dot);

using FloatTypes = ::testing::Types<std::pair<double, double>, std::pair<double, float>>;
INSTANTIATE_TYPED_TEST_SUITE_P(Float, TestArrayHandlerSparse, FloatTypes);

TEST(ArrayHandlerSparse, copy) {
  using X = std::map<size_t, double>;
  using Y = std::map<size_t, float>;
  ArrayHandlerSparse<X, Y> handler{};
  auto y = Y{{0, 1}, {4, 2}, {5, 3}, {11, 4}};
  auto x = handler.copy(y);
  ASSERT_EQ(x.size(), y.size());
  auto iy = y.begin();
  auto ix = x.begin();
  for (auto i = 0; i < x.size(); ++i, ++iy, ++ix) {
    EXPECT_DOUBLE_EQ(iy->second, ix->second) << "i = " << std::to_string(i);
    EXPECT_EQ(iy->first, ix->first) << "i = " << std::to_string(i);
  }
}
