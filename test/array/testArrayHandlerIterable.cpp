#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <deque>
#include <molpro/linalg/array/ArrayHandlerIterable.h>

using molpro::linalg::array::ArrayHandlerIterable;

using ::testing::ContainerEq;
using ::testing::Each;

TEST(ArrayHandlerIterable, constructor) { ArrayHandlerIterable<std::vector<int>> handler{}; }

TEST(ArrayHandlerIterable, lazy_dot) {
  using value_type = int;
  static const int N = 2;
  static const int dim = 5;
  ArrayHandlerIterable<std::vector<value_type>> handler{};
  auto xx = std::vector<std::vector<value_type>>(N);
  auto yy = xx;
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

TEST(ArrayHandlerIterable, lazy_axpy) {
  using value_type = int;
  static const int N = 2;
  static const int dim = 5;
  static const value_type alpha = 3;
  static const value_type xval = 2;
  static const value_type yval = 5;
  ArrayHandlerIterable<std::vector<value_type>> handler{};
  auto xx = std::vector<std::vector<value_type>>(N, std::vector<value_type>(dim, xval));
  auto yy = std::vector<std::vector<value_type>>(N, std::vector<value_type>(dim, yval));
  auto ref_axpy = std::vector<std::vector<value_type>>(N, std::vector<value_type>(dim, xval + alpha * yval));
  {
    auto h = handler.lazy_handle();
    for (size_t i = 0; i < N; ++i) {
      h.axpy(xx[i], alpha, yy[i]);
    }
    EXPECT_FALSE(h.invalid());
    for (const auto& x : xx)
      EXPECT_THAT(x, Each(xval));
  }
  for (size_t i = 0; i < xx.size(); ++i)
    EXPECT_THAT(xx[i], ContainerEq(ref_axpy[i]));
}

TEST(ArrayHandlerIterable, copy_same) {
  using X = std::vector<double>;
  using Y = X;
  ArrayHandlerIterable<X, Y> handler{};
  auto x = X{1, 2, 3, 4};
  auto y = handler.copy(x);
}

TEST(ArrayHandlerIterable, copyL) {
  using X = std::vector<double>;
  using Y = std::deque<int>;
  ArrayHandlerIterable<X, Y> handler{};
  auto x = X{1, 2, 3, 4};
  auto y = handler.copy(x);
}

TEST(ArrayHandlerIterable, copyR) {
  using X = std::vector<double>;
  using Y = std::deque<int>;
  ArrayHandlerIterable<X, Y> handler{};
  auto y = Y{1, 2, 3, 4};
  auto x = handler.copy(y);
}
