#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/array/ArrayHandler.h>
#include <numeric>
#include <random>

using molpro::linalg::array::util::OperationRegister;
using molpro::linalg::array::util::RefEqual;
using molpro::linalg::array::util::remove_duplicates;

using ::testing::ContainerEq;

// template this based on number of arrays
class ArrayHandlerInputF : public ::testing::Test {
public:
  using value_type = long int;
  using X = std::vector<value_type>;
  using Y = std::vector<int>;
  ArrayHandlerInputF() : x_arrays(na), y_arrays(na), scalars(na * na) {}

  std::vector<X> x_arrays;
  std::vector<Y> y_arrays;
  std::vector<value_type> scalars;
  static const size_t na = 3; // number of arrays
};

struct OperationRegisterInput : public ::testing::Test {
  using X = int;
  using Y = double;
  OperationRegisterInput() : x(size), y(size) {
    std::iota(begin(x), end(x), 0);
    std::iota(begin(x), end(x), 1);
  }

  static constexpr int size = 5;
  std::vector<X> x;
  std::vector<Y> y;
};

TEST_F(OperationRegisterInput, no_priority) {
  auto op_register = OperationRegister<X, Y>{};
  decltype(op_register.m_register) register_ref{};
  for (size_t i = 0; i < size; ++i)
    for (size_t j = 0; j < size; ++j) {
      op_register.push(x[i], y[i]);
      register_ref.emplace_back(x[i], y[i]);
    }
  EXPECT_THAT(op_register.m_register, ContainerEq(register_ref));
}

TEST_F(OperationRegisterInput, prioritize_first_element) {
  auto op_register = OperationRegister<X, Y>{};
  decltype(op_register.m_register) register_ref{};
  for (size_t i = 0; i < size; ++i)
    for (size_t j = 0; j < size; ++j) {
      op_register.push<0, std::equal_to<X>>(x[i], y[j], {});
      register_ref.emplace_back(x[i], y[i]);
    }
  EXPECT_THAT(op_register.m_register, ContainerEq(register_ref));
}

TEST_F(OperationRegisterInput, prioritize_second_element) {
  auto op_register = OperationRegister<X, Y>{};
  decltype(op_register.m_register) register_ref{};
  for (size_t i = 0; i < size; ++i)
    for (size_t j = 0; j < size; ++j) {
      op_register.push<1, std::equal_to<Y>>(x[j], y[i], {});
    }
  for (size_t i = 0; i < size; ++i)
    for (size_t j = 0; j < size; ++j) {
      register_ref.emplace_back(x[j], y[i]);
    }
  EXPECT_THAT(op_register.m_register, ContainerEq(register_ref));
}

TEST(ArrayHandlerUtil, remove_duplicates_small) {
  using X = double;
  using Y = int;
  using S = size_t;
  std::vector<X> x{0, 1};
  std::vector<Y> y{0, 1};
  std::vector<S> s{0, 1, 2, 3, 4};
  std::list<std::tuple<std::reference_wrapper<X>, std::reference_wrapper<Y>, std::reference_wrapper<S>>> op_register{
      {x[0], y[0], s[0]}, {x[0], y[1], s[1]}, {x[1], y[1], s[2]},
      {x[0], y[0], s[3]}, {x[1], y[0], s[4]}, {x[1], y[0], s[4]}};
  std::vector<std::reference_wrapper<X>> ref_unique_x{x[0], x[1]};
  std::vector<std::reference_wrapper<Y>> ref_unique_y{y[0], y[1]};
  std::vector<std::reference_wrapper<S>> ref_unique_s{s[0], s[1], s[2], s[3], s[4]};
  std::vector<std::tuple<size_t, size_t, size_t>> ref_op_register_ind{{0, 0, 0}, {0, 1, 1}, {1, 1, 2},
                                                                      {0, 0, 3}, {1, 0, 4}, {1, 0, 4}};
  std::vector<void*> ref_addr_x;
  std::vector<void*> ref_addr_y;
  std::vector<void*> ref_addr_s;
  auto addr_of = [](auto& el) { return std::addressof(el.get()); };
  std::transform(begin(ref_unique_x), end(ref_unique_x), std::back_inserter(ref_addr_x), addr_of);
  std::transform(begin(ref_unique_y), end(ref_unique_y), std::back_inserter(ref_addr_y), addr_of);
  std::transform(begin(ref_unique_s), end(ref_unique_s), std::back_inserter(ref_addr_s), addr_of);
  auto res = remove_duplicates<std::reference_wrapper<X>, std::reference_wrapper<Y>, std::reference_wrapper<S>,
                               RefEqual<X>, RefEqual<Y>, RefEqual<S>>(op_register, {}, {}, {});
  auto unique_x = std::get<1>(res);
  auto unique_y = std::get<2>(res);
  auto unique_s = std::get<3>(res);
  std::vector<void*> addr_x;
  std::vector<void*> addr_y;
  std::vector<void*> addr_s;
  std::transform(begin(unique_x), end(unique_x), std::back_inserter(addr_x), addr_of);
  std::transform(begin(unique_y), end(unique_y), std::back_inserter(addr_y), addr_of);
  std::transform(begin(unique_s), end(unique_s), std::back_inserter(addr_s), addr_of);
  EXPECT_THAT(std::get<0>(res), ContainerEq(ref_op_register_ind));
  EXPECT_THAT(addr_x, ContainerEq(ref_addr_x));
  EXPECT_THAT(addr_y, ContainerEq(ref_addr_y));
  EXPECT_THAT(addr_s, ContainerEq(ref_addr_s));
}
