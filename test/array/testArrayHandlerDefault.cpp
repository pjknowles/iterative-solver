#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <numeric>
#include <utility>

#include <molpro/linalg/array/ArrayHandlerDefault.h>

using molpro::linalg::array::ArrayHandlerDefault;

using ::testing::ContainerEq;
using ::testing::Each;

struct DummyArray {
  using value_type = double;
};

struct DefaultArray {
  using value_type = double;
  DefaultArray() = default;
  DefaultArray(const DummyArray &) : DefaultArray() { ++count_copy; };
  void scal(value_type alpha) { ++count_scal; }
  void fill(value_type alpha) { ++count_fill; }
  void axpy(value_type alpha, const DummyArray &x) { ++count_axpy; }
  value_type dot(const DummyArray &y) const {
    ++count_dot;
    return 0.;
  }
  std::map<size_t, value_type> select_max_dot(size_t n, const DummyArray &x) const {
    ++count_select_max_dot;
    return {};
  }
  int count_copy = 0;
  int count_scal = 0;
  int count_fill = 0;
  int count_axpy = 0;
  mutable int count_dot = 0;
  mutable int count_select_max_dot = 0;
};

TEST(TestArrayHandlerDefault, constructor) {
  auto a = DummyArray{};
  auto b = DefaultArray{};
  auto handler = ArrayHandlerDefault<DefaultArray, DummyArray>();
  auto c = handler.copy(a);
  EXPECT_EQ(c.count_copy, 1);
  handler.scal(1., b);
  EXPECT_EQ(b.count_scal, 1);
  handler.fill(1., b);
  EXPECT_EQ(b.count_fill, 1);
  handler.axpy(1., a, b);
  EXPECT_EQ(b.count_axpy, 1);
  handler.dot(b, a);
  EXPECT_EQ(b.count_dot, 1);
  handler.select_max_dot(10, b, a);
  EXPECT_EQ(b.count_select_max_dot, 1);
}