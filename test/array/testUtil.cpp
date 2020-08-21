#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/array/util/select_max_dot.h>

#include <complex>

using molpro::linalg::array::util::select_max_dot;
using ::testing::ContainerEq;

template <class T>
struct Test_select_max_dot : ::testing::Test {};

TYPED_TEST_SUITE_P(Test_select_max_dot);

TYPED_TEST_P(Test_select_max_dot, std_vector) {
  using std::abs;
  using value_type_L = typename TypeParam::first_type;
  using value_type_R = typename TypeParam::second_type;
  using value_type = decltype(value_type_L{} * value_type_R{});
  using value_type_abs = decltype(abs(value_type{}));
  auto x = std::vector<value_type_L>{1, -2, 1, 0, 3, 0, -4, 1};
  auto y = std::vector<value_type_R>{1, 1, 1, 1, 1, 1, 1, 1};
  auto ref_result = std::map<size_t, value_type_abs>{{6, 4}, {4, 3}, {1, 2}};
  auto select = select_max_dot<decltype(x), decltype(y), value_type, value_type_abs>(ref_result.size(), x, y);
  ASSERT_THAT(select, ContainerEq(ref_result));
}

REGISTER_TYPED_TEST_SUITE_P(Test_select_max_dot, std_vector);

using IntTypes = ::testing::Types<std::pair<int, int>, std::pair<int, short>, std::pair<short, int>>;
using FloatTypes = ::testing::Types<std::pair<double, double>, std::pair<double, float>, std::pair<float, double>>;
using MixTypes = ::testing::Types<std::pair<double, int>, std::pair<int, double>>;
using ComplexTypes = ::testing::Types<std::pair<std::complex<double>, std::complex<double>>>;
INSTANTIATE_TYPED_TEST_SUITE_P(Int, Test_select_max_dot, IntTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(Float, Test_select_max_dot, FloatTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(Mix, Test_select_max_dot, MixTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(Complex, Test_select_max_dot, ComplexTypes);
