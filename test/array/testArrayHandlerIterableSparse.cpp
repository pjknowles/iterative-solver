#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/array/ArrayHandlerIterableSparse.h>

using molpro::linalg::array::ArrayHandlerIterableSparse;

using ::testing::DoubleEq;
using ::testing::Pointwise;

template <class T>
struct TestArrayHandlerIterableSparse : ::testing::Test {};

TYPED_TEST_SUITE_P(TestArrayHandlerIterableSparse);

TYPED_TEST_P(TestArrayHandlerIterableSparse, constructor) {
  auto handler = ArrayHandlerIterableSparse<std::vector<typename TypeParam::first_type>,
                                            std::map<size_t, typename TypeParam::second_type>>{};
}

TYPED_TEST_P(TestArrayHandlerIterableSparse, dot) {
  auto x = std::vector<typename TypeParam::first_type>{1, 1.0, 1, 2.0, 1, 3.0, 0, 1, 1, 1, 1, 4.0};
  auto y = std::map<size_t, typename TypeParam::second_type>{{1, 1.0}, {3, 2.0}, {6, 3.0}, {11, 4.0}};
  auto handler = ArrayHandlerIterableSparse<decltype(x), decltype(y)>{};
  auto ref = typename decltype(handler)::value_type(1. + 2. * 2. + 4. * 4.);
  auto result = handler.dot(x, y);
  ASSERT_DOUBLE_EQ(result, ref);
}

REGISTER_TYPED_TEST_SUITE_P(TestArrayHandlerIterableSparse, constructor, dot);

using FloatTypes = ::testing::Types<std::pair<double, double>, std::pair<double, float>>;
INSTANTIATE_TYPED_TEST_SUITE_P(Float, TestArrayHandlerIterableSparse, FloatTypes);
