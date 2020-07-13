#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/array/ArrayHandler.h>
#include <numeric>
#include <random>

using molpro::linalg::array::util::RegisterOperation;

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

template <class OP_Type> void test_op_registers(const std::list<OP_Type>& reg1, const std::list<OP_Type>& reg2) {
  auto ops_equal = [](const OP_Type& op1, const OP_Type& op2) {
    return std::addressof(std::get<0>(op1).get()) == std::addressof(std::get<0>(op2).get()) &&
           std::addressof(std::get<1>(op1).get()) == std::addressof(std::get<1>(op2).get()) &&
           std::addressof(std::get<2>(op1).get()) == std::addressof(std::get<2>(op2).get());
  };
  ASSERT_EQ(reg1.size(), reg2.size());
  for (size_t i = 0; i < reg1.size(); ++i) {
    EXPECT_TRUE(ops_equal(*std::next(reg1.cbegin(), i), *std::next(reg2.cbegin(), i))) << " i = " << std::to_string(i);
  }
}

struct RegisterOperationLeft : ArrayHandlerInputF {
  using OP_Type =
      std::tuple<std::reference_wrapper<X>, std::reference_wrapper<const Y>, std::reference_wrapper<value_type>>;
  RegisterOperationLeft() {
    for (size_t i = 0, ij = 0; i < x_arrays.size(); ++i)
      for (size_t j = 0; j < y_arrays.size(); ++j, ++ij)
        op_register_ref.emplace_back(x_arrays[i], y_arrays[j], scalars[ij]);
  }

  std::list<OP_Type> op_register_ref{};
};

TEST_F(RegisterOperationLeft, no_priority) {
  std::list<OP_Type> op_register{};
  for (size_t i = 0, ij = 0; i < x_arrays.size(); ++i)
    for (size_t j = 0; j < y_arrays.size(); ++j, ++ij)
      RegisterOperation<OP_Type, -1>()(op_register, {x_arrays[i], y_arrays[j], scalars[ij]});
  test_op_registers(op_register, op_register_ref);
}

TEST_F(RegisterOperationLeft, no_reorder) {
  std::list<OP_Type> op_register{};
  for (size_t i = 0, ij = 0; i < x_arrays.size(); ++i)
    for (size_t j = 0; j < y_arrays.size(); ++j, ++ij)
      RegisterOperation<OP_Type, 0>()(op_register, {x_arrays[i], y_arrays[j], scalars[ij]});
  test_op_registers(op_register, op_register_ref);
}

TEST_F(RegisterOperationLeft, reverse_loop_order) {
  std::list<OP_Type> op_register{};
  for (size_t j = 0; j < y_arrays.size(); ++j)
    for (size_t i = 0; i < x_arrays.size(); ++i)
      RegisterOperation<OP_Type, 0>()(op_register, {x_arrays[i], y_arrays[j], scalars[i * x_arrays.size() + j]});
  test_op_registers(op_register, op_register_ref);
}

TEST_F(RegisterOperationLeft, random_order) {
  std::list<OP_Type> op_register{};
  std::vector<std::pair<size_t, size_t>> rand_ind{{1, 2}, {0, 1}, {0, 0}, {2, 0}, {0, 2},
                                                  {1, 0}, {1, 1}, {2, 1}, {2, 2}};
  std::vector<std::pair<size_t, size_t>> rand_ind_ref{{1, 2}, {1, 0}, {1, 1}, {0, 1}, {0, 0},
                                                      {0, 2}, {2, 0}, {2, 1}, {2, 2}};
  op_register_ref.clear();
  for (const auto& ind : rand_ind_ref) {
    size_t i, j;
    std::tie(i, j) = ind;
    op_register_ref.emplace_back(x_arrays[i], y_arrays[j], scalars[i * y_arrays.size() + j]);
  }
  for (const auto& ind : rand_ind) {
    size_t i, j;
    std::tie(i, j) = ind;
    RegisterOperation<OP_Type, 0>()(op_register, {x_arrays[i], y_arrays[j], scalars[i * y_arrays.size() + j]});
  }
  test_op_registers(op_register, op_register_ref);
}

struct RegisterOperationRight : ArrayHandlerInputF {
  using OP_Type =
      std::tuple<std::reference_wrapper<X>, std::reference_wrapper<const Y>, std::reference_wrapper<value_type>>;
  RegisterOperationRight() {
    for (size_t i = 0, ij = 0; i < y_arrays.size(); ++i)
      for (size_t j = 0; j < x_arrays.size(); ++j, ++ij)
        op_register_ref.emplace_back(x_arrays[j], y_arrays[i], scalars[ij]);
  }

  std::list<OP_Type> op_register_ref{};
};

TEST_F(RegisterOperationRight, no_priority) {
  std::list<OP_Type> op_register{};
  for (size_t i = 0, ij = 0; i < y_arrays.size(); ++i)
    for (size_t j = 0; j < x_arrays.size(); ++j, ++ij)
      RegisterOperation<OP_Type, -1>()(op_register, {x_arrays[j], y_arrays[i], scalars[ij]});
  test_op_registers(op_register, op_register_ref);
}

TEST_F(RegisterOperationRight, no_reorder) {
  std::list<OP_Type> op_register{};
  for (size_t i = 0, ij = 0; i < y_arrays.size(); ++i)
    for (size_t j = 0; j < x_arrays.size(); ++j, ++ij)
      RegisterOperation<OP_Type, 1>()(op_register, {x_arrays[j], y_arrays[i], scalars[ij]});
  test_op_registers(op_register, op_register_ref);
}

TEST_F(RegisterOperationRight, reverse_loop_order) {
  std::list<OP_Type> op_register{};
  for (size_t i = 0; i < x_arrays.size(); ++i)
    for (size_t j = 0; j < y_arrays.size(); ++j)
      RegisterOperation<OP_Type, 1>()(op_register, {x_arrays[i], y_arrays[j], scalars[j * x_arrays.size() + i]});
  test_op_registers(op_register, op_register_ref);
}

TEST_F(RegisterOperationRight, random_order) {
  std::list<OP_Type> op_register{};
  std::vector<std::pair<size_t, size_t>> rand_ind{{1, 2}, {0, 1}, {0, 0}, {2, 0}, {0, 2},
                                                  {1, 0}, {1, 1}, {2, 1}, {2, 2}};
  std::vector<std::pair<size_t, size_t>> rand_ind_ref{{1, 2}, {0, 2}, {2, 2}, {0, 1}, {1, 1},
                                                      {2, 1}, {0, 0}, {2, 0}, {1, 0}};
  op_register_ref.clear();
  for (const auto& ind : rand_ind_ref) {
    size_t i, j;
    std::tie(i, j) = ind;
    op_register_ref.emplace_back(x_arrays[i], y_arrays[j], scalars[i * y_arrays.size() + j]);
  }
  for (const auto& ind : rand_ind) {
    size_t i, j;
    std::tie(i, j) = ind;
    RegisterOperation<OP_Type, 1>()(op_register, {x_arrays[i], y_arrays[j], scalars[i * y_arrays.size() + j]});
  }
  test_op_registers(op_register, op_register_ref);
}

// TODO test remove_duplicates
