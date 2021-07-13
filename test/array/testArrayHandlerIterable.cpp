#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <complex>
#include <deque>
#include <molpro/linalg/array/ArrayHandlerIterable.h>

using molpro::linalg::array::ArrayHandlerIterable;

using ::testing::ContainerEq;
using ::testing::DoubleEq;
using ::testing::Each;
using ::testing::Pointwise;

template <class T>
struct TestArrayHandlerIterable : ::testing::Test {};

TYPED_TEST_SUITE_P(TestArrayHandlerIterable);

TYPED_TEST_P(TestArrayHandlerIterable, constructor) {
  ArrayHandlerIterable<std::vector<typename TypeParam::first_type>, std::vector<typename TypeParam::second_type>>
      handler{};
}

TYPED_TEST_P(TestArrayHandlerIterable, lazy_dot) {
  ArrayHandlerIterable<std::vector<typename TypeParam::first_type>, std::vector<typename TypeParam::second_type>>
      handler{};
  using value_type_L = typename decltype(handler)::value_type_L;
  using value_type_R = typename decltype(handler)::value_type_R;
  using value_type = typename decltype(handler)::value_type;
  static const int N = 3;
  static const int dim = 5;
  auto xx = std::vector<std::vector<value_type_L>>(N);
  auto yy = std::vector<std::vector<value_type_R>>(N);
  auto result = std::vector<value_type>(N * N);
  auto ref_dot = std::vector<value_type>(N * N);
  for (size_t i = 0; i < N; ++i) {
    xx[i].resize(dim);
    yy[i].resize(dim);
    std::iota(begin(xx[i]), end(xx[i]), i);
    std::iota(begin(yy[i]), end(yy[i]), 10 * i);
  }
  {
    auto h = handler.lazy_handle();
    for (size_t i = 0, ij = 0; i < N; ++i) {
      for (size_t j = 0; j < N; ++j, ++ij) {
        ref_dot[ij] = handler.dot(xx[i], yy[j]);
        h.dot(xx[i], yy[j], result[ij]);
      }
    }
    EXPECT_FALSE(h.invalid());
    EXPECT_THAT(result, Each(0));
  }
  EXPECT_THAT(result, ContainerEq(ref_dot));
}

TYPED_TEST_P(TestArrayHandlerIterable, select_max_dot) {
  ArrayHandlerIterable<std::vector<typename TypeParam::first_type>, std::vector<typename TypeParam::second_type>>
      handler{};
  using value_type_L = typename decltype(handler)::value_type_L;
  using value_type_R = typename decltype(handler)::value_type_R;
//  using value_type = typename decltype(handler)::value_type;
  using value_type_abs = typename decltype(handler)::value_type_abs;
  auto x = std::vector<value_type_L>{1, -2, 1, 0, 3, 0, -4, 1};
  auto y = std::vector<value_type_R>{1, 1, 1, 1, 1, 1, 1, 1};
  auto ref_result = std::map<size_t, value_type_abs>{{6, 4}, {4, 3}, {1, 2}};
  auto select = handler.select_max_dot(ref_result.size(), x, y);
  ASSERT_THAT(select, ContainerEq(ref_result));
}

TYPED_TEST_P(TestArrayHandlerIterable, lazy_axpy) {
  ArrayHandlerIterable<std::vector<typename TypeParam::first_type>, std::vector<typename TypeParam::second_type>>
      handler{};
  using value_type_L = typename decltype(handler)::value_type_L;
  using value_type_R = typename decltype(handler)::value_type_R;
  using value_type = typename decltype(handler)::value_type;
  static const int N = 3;
  static const int dim = 5;
  static const value_type alpha = 3;
  static const value_type_R xval = 2;
  static const value_type_L yval = 5;
  auto xx = std::vector<std::vector<value_type_R>>(N, std::vector<value_type_R>(dim, xval));
  auto yy = std::vector<std::vector<value_type_L>>(N, std::vector<value_type_L>(dim, yval));
  auto ref_axpy = std::vector<std::vector<value_type_L>>(N, std::vector<value_type_L>(dim, alpha * xval + yval));
  {
    auto h = handler.lazy_handle();
    for (size_t i = 0; i < N; ++i) {
      h.axpy(alpha, xx[i], yy[i]);
    }
    EXPECT_FALSE(h.invalid());
    for (const auto& y : yy)
      EXPECT_THAT(y, Each(yval));
  }
  for (size_t i = 0; i < xx.size(); ++i)
    EXPECT_THAT(yy[i], ContainerEq(ref_axpy[i]));
}

TYPED_TEST_P(TestArrayHandlerIterable, lazy_axpy_lazy_off) {
  ArrayHandlerIterable<std::vector<typename TypeParam::first_type>, std::vector<typename TypeParam::second_type>>
      handler{};
  using value_type_L = typename decltype(handler)::value_type_L;
  using value_type_R = typename decltype(handler)::value_type_R;
  using value_type = typename decltype(handler)::value_type;
  static const int N = 3;
  static const int dim = 5;
  static const value_type alpha = 3;
  static const value_type_R xval = 2;
  static const value_type_L yval = 5;
  auto xx = std::vector<std::vector<value_type_R>>(N, std::vector<value_type_R>(dim, xval));
  auto yy = std::vector<std::vector<value_type_L>>(N, std::vector<value_type_L>(dim, yval));
  auto ref_axpy = std::vector<std::vector<value_type_L>>(N, std::vector<value_type_L>(dim, alpha * xval + yval));
  {
    auto h = handler.lazy_handle();
    EXPECT_FALSE(h.is_off());
    h.off();
    for (size_t i = 0; i < N; ++i) {
      h.axpy(alpha, xx[i], yy[i]);
    }
    EXPECT_FALSE(h.invalid());
    EXPECT_TRUE(h.is_off());
    for (size_t i = 0; i < xx.size(); ++i)
      EXPECT_THAT(yy[i], ContainerEq(ref_axpy[i]));
  }
}

REGISTER_TYPED_TEST_SUITE_P(TestArrayHandlerIterable, constructor, lazy_dot, select_max_dot, lazy_axpy,
                            lazy_axpy_lazy_off);

using IntTypes = ::testing::Types<std::pair<int, int>, std::pair<int, short>, std::pair<short, int>>;
using FloatTypes = ::testing::Types<std::pair<double, double>, std::pair<double, float>, std::pair<float, double>>;
using MixTypes = ::testing::Types<std::pair<double, int>, std::pair<int, double>>;
using ComplexTypes = ::testing::Types<std::pair<std::complex<double>, std::complex<double>>>;
INSTANTIATE_TYPED_TEST_SUITE_P(Int, TestArrayHandlerIterable, IntTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(Float, TestArrayHandlerIterable, FloatTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(Mix, TestArrayHandlerIterable, MixTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(Complex, TestArrayHandlerIterable, ComplexTypes);

TEST(ArrayHandlerIterable, copy_same) {
  using X = std::vector<double>;
  using Y = X;
  ArrayHandlerIterable<X, Y> handler{};
  auto y = Y{1, 2, 3, 4};
  auto x = handler.copy(y);
  EXPECT_THAT(x, Pointwise(DoubleEq(), y));
}

TEST(ArrayHandlerIterable, copyL) {
  using X = std::vector<double>;
  using Y = std::deque<int>;
  ArrayHandlerIterable<X, Y> handler{};
  auto y = Y{1, 2, 3, 4};
  auto x = handler.copy(y);
  EXPECT_THAT(x, Pointwise(DoubleEq(), y));
}

TEST(ArrayHandlerIterable, copyR) {
  using X = std::vector<double>;
  using Y = std::deque<int>;
  ArrayHandlerIterable<X, Y> handler{};
  auto y = Y{1, 2, 3, 4};
  auto x = handler.copy(y);
  EXPECT_THAT(x, Pointwise(DoubleEq(), y));
}
