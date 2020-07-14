#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/array/ArrayHandlerFactory.h>

#include <utility>
using molpro::linalg::array::ArrayHandler;
using molpro::linalg::array::ArrayHandlerFactory;
using molpro::linalg::array::ArrayHandlerIterable;

TEST(ArrayHandlerFactory, iterable) {
  using X = std::vector<int>;
  using Y = std::vector<long int>;
  auto handle = ArrayHandlerFactory<X, Y>::create();
  using handleType = std::remove_reference_t<decltype(*handle)>;
  constexpr bool typeHandler = std::is_base_of<ArrayHandler<X, Y>, handleType>::value;
  constexpr bool typeIterable = std::is_same<ArrayHandlerIterable<X, Y>, handleType>::value;
  EXPECT_TRUE(typeHandler);
  EXPECT_TRUE(typeIterable);
}

TEST(ArrayHandlerFactory, DummyIterativeSolver) {
  using X = std::vector<int>;
  using Y = std::vector<long int>;
  struct IterativeSolver {
    IterativeSolver(std::shared_ptr<ArrayHandler<X, Y>> handler = ArrayHandlerFactory<X, Y>::create())
        : m_handler{std::move(handler)} {}
    std::shared_ptr<ArrayHandler<X, Y>> m_handler;
  };
  auto solver = IterativeSolver();
}
