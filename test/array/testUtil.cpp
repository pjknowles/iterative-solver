#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/array/util/select_max_dot.h>

#include <complex>

using molpro::linalg::array::util::select_max_dot;
using molpro::linalg::array::util::select_max_dot_iter_sparse;
using molpro::linalg::array::util::select_max_dot_sparse;
using ::testing::ContainerEq;

template <class T>
struct Test_select_max_dot : ::testing::Test {};

TYPED_TEST_SUITE_P(Test_select_max_dot);

TYPED_TEST_P(Test_select_max_dot, iterable) {
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

TYPED_TEST_P(Test_select_max_dot, iterable_sparse) {
  using std::abs;
  using value_type_L = typename TypeParam::first_type;
  using value_type_R = typename TypeParam::second_type;
  using value_type = decltype(value_type_L{} * value_type_R{});
  using value_type_abs = decltype(abs(value_type{}));
  auto x = std::vector<value_type_L>{1, -2, 1, 0, 3, 0, -4, 1};
  auto y = std::map<size_t, value_type_R>{{0, 1}, {1, 1}, {2, 1}, {4, 1}, {6, 1}};
  auto ref_result = std::map<size_t, value_type_abs>{{6, 4}, {4, 3}, {1, 2}};
  auto select =
      select_max_dot_iter_sparse<decltype(x), decltype(y), value_type, value_type_abs>(ref_result.size(), x, y);
  ASSERT_THAT(select, ContainerEq(ref_result));
}

TYPED_TEST_P(Test_select_max_dot, sparse) {
  using std::abs;
  using value_type_L = typename TypeParam::first_type;
  using value_type_R = typename TypeParam::second_type;
  using value_type = decltype(value_type_L{} * value_type_R{});
  using value_type_abs = decltype(abs(value_type{}));
  auto x = std::map<size_t, value_type_L>{{1, -2}, {3, 1}, {4, 3}, {6, -4}};
  auto y = std::map<size_t, value_type_R>{{0, 1}, {1, 1}, {2, 1}, {4, 1}, {6, 1}};
  auto ref_result = std::map<size_t, value_type_abs>{{6, 4}, {4, 3}, {1, 2}};
  auto select = select_max_dot_sparse<decltype(x), decltype(y), value_type, value_type_abs>(ref_result.size(), x, y);
  ASSERT_THAT(select, ContainerEq(ref_result));
  select = select_max_dot_sparse<decltype(x), decltype(y), value_type, value_type_abs>(ref_result.size() + 1, x, y);
  ASSERT_THAT(select, ContainerEq(ref_result)) << "sparse should be able to return a smaller selection than asked for";
}

REGISTER_TYPED_TEST_SUITE_P(Test_select_max_dot, iterable, iterable_sparse, sparse);

using Types = ::testing::Types<std::pair<int, int>, std::pair<double, double>,
                               std::pair<std::complex<double>, std::complex<double>>,
                               std::pair<double, std::complex<double>>, std::pair<std::complex<double>, double>>;
INSTANTIATE_TYPED_TEST_SUITE_P(Types, Test_select_max_dot, Types);
