#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/array/ArrayHandlerIterable.h>

using molpro::linalg::array::ArrayHandlerIterable;

using ::testing::ContainerEq;

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
  }
  EXPECT_THAT(result, ContainerEq(ref_dot));
}

//TODO test axpy

