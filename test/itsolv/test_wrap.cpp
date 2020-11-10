#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/itsolv/wrap.h>
#include <molpro/linalg/itsolv/wrap_util.h>

#include <list>

using molpro::linalg::itsolv::cwrap;
using molpro::linalg::itsolv::find_ref;
using molpro::linalg::itsolv::remove_elements;
using molpro::linalg::itsolv::wrap;
using ::testing::Eq;
using ::testing::Pointwise;

TEST(wrap_util, find_ref_empty) {
  auto params = std::vector<int>{};
  auto result = find_ref(wrap(params), begin(params), end(params));
  ASSERT_TRUE(result.empty());
}

TEST(wrap_util, find_ref_all) {
  auto params = std::vector<int>{1, 2, 3, 4};
  auto indices = find_ref(wrap(params), begin(params), end(params));
  auto ref_indices = std::vector<size_t>{0, 1, 2, 3};
  ASSERT_FALSE(indices.empty());
  ASSERT_THAT(indices, Pointwise(Eq(), ref_indices));
}

TEST(wrap_util, find_ref_all_const) {
  const auto params = std::vector<int>{1, 2, 3, 4};
  auto indices = find_ref(wrap(params), begin(params), end(params));
  auto ref_indices = std::vector<size_t>{0, 1, 2, 3};
  ASSERT_FALSE(indices.empty());
  ASSERT_THAT(indices, Pointwise(Eq(), ref_indices));
}

TEST(wrap_util, find_ref_some) {
  auto params = std::vector<int>{1, 2, 3, 4, 5};
  auto wparams = wrap(params);
  wparams.erase(wparams.begin());
  wparams.erase(wparams.begin() + 1);
  wparams.erase(wparams.begin() + 2);
  auto indices = find_ref(wparams, begin(params), end(params));
  auto ref_indices = std::vector<size_t>{1, 3};
  ASSERT_FALSE(indices.empty());
  ASSERT_THAT(indices, Pointwise(Eq(), ref_indices));
}

TEST(wrap_util, remove_elements) {
  auto params = std::vector<int>{1, 2, 3, 4, 5, 6};
  auto indices = std::vector<size_t>{0, 2, 3, 5};
  auto params_ref = std::vector<int>{2, 5};
  auto result = remove_elements(params, indices);
  ASSERT_EQ(result.size(), params_ref.size());
  ASSERT_THAT(result, Pointwise(Eq(), params_ref));
}

TEST(wrap, iterable_container__list) {
  auto params = std::list<int>{1, 2, 3, 4, 5, 6};
  auto wparams = wrap<int>(params);
  ASSERT_FALSE(std::is_const_v<decltype(wparams)::value_type::type>);
  ASSERT_EQ(wparams.size(), params.size());
  ASSERT_THAT(wparams, Pointwise(Eq(), params));
}

TEST(wrap, iterable_container__const_list) {
  const auto params = std::list<int>{1, 2, 3, 4, 5, 6};
  auto wparams = wrap<int>(params);
  ASSERT_TRUE(std::is_const_v<decltype(wparams)::value_type::type>);
  ASSERT_EQ(wparams.size(), params.size());
  ASSERT_THAT(wparams, Pointwise(Eq(), params));
}

TEST(cwrap, iterable_container__list) {
  const auto params = std::list<int>{1, 2, 3, 4, 5, 6};
  auto wparams = cwrap<int>(params);
  ASSERT_TRUE(std::is_const_v<decltype(wparams)::value_type::type>);
  ASSERT_EQ(wparams.size(), params.size());
  ASSERT_THAT(wparams, Pointwise(Eq(), params));
}

TEST(cwrap, iterable_container__const_list) {
  const auto params = std::list<int>{1, 2, 3, 4, 5, 6};
  auto wparams = cwrap<int>(params);
  ASSERT_TRUE(std::is_const_v<decltype(wparams)::value_type::type>);
  ASSERT_EQ(wparams.size(), params.size());
  ASSERT_THAT(wparams, Pointwise(Eq(), params));
}

TEST(cwrap, VecRef_to_CVecRef) {
  auto params = std::vector<int>{1, 2, 3, 4, 5, 6};
  auto wparams = wrap(params);
  auto cwparams = cwrap(wparams);
  ASSERT_TRUE(std::is_const_v<decltype(cwparams)::value_type::type>);
  ASSERT_THAT(wparams, Pointwise(Eq(), params));
  ASSERT_THAT(cwparams, Pointwise(Eq(), params));
  ASSERT_THAT(cwparams, Pointwise(Eq(), wparams));
}

TEST(cwrap, const_VecRef_to_CVecRef) {
  auto params = std::vector<int>{1, 2, 3, 4, 5, 6};
  const auto wparams = wrap(params);
  auto cwparams = cwrap(wparams);
  ASSERT_TRUE(std::is_const_v<decltype(cwparams)::value_type::type>);
  ASSERT_THAT(wparams, Pointwise(Eq(), params));
  ASSERT_THAT(cwparams, Pointwise(Eq(), params));
  ASSERT_THAT(cwparams, Pointwise(Eq(), wparams));
}
