#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/itsolv/wrap.h>

using molpro::linalg::itsolv::find_ref;
using molpro::linalg::itsolv::wrap;
using ::testing::Eq;
using ::testing::Pointwise;

TEST(wrap, find_ref_empty) {
  auto params = std::vector<int>{};
  auto result = find_ref(wrap(params), begin(params), end(params));
  ASSERT_TRUE(result.empty());
}

TEST(wrap, find_ref_all) {
  auto params = std::vector<int>{1, 2, 3, 4};
  auto indices = find_ref(wrap(params), begin(params), end(params));
  auto ref_indices = std::vector<size_t>{0, 1, 2, 3};
  ASSERT_FALSE(indices.empty());
  ASSERT_THAT(indices, Pointwise(Eq(), ref_indices));
}

TEST(wrap, find_ref_all_const) {
  const auto params = std::vector<int>{1, 2, 3, 4};
  auto indices = find_ref(wrap(params), begin(params), end(params));
  auto ref_indices = std::vector<size_t>{0, 1, 2, 3};
  ASSERT_FALSE(indices.empty());
  ASSERT_THAT(indices, Pointwise(Eq(), ref_indices));
}

TEST(wrap, find_ref_some) {
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
