#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/itsolv/util.h>

using molpro::linalg::itsolv::util::construct_zeroed_copy;
using molpro::linalg::itsolv::util::delete_parameters;
using molpro::linalg::itsolv::util::is_iota;
using molpro::linalg::itsolv::util::StringFacet;

TEST(itsolv_util, is_iota_null) {
  auto vec = std::vector<int>{};
  ASSERT_TRUE(is_iota(vec.begin(), vec.end(), 0));
}

TEST(itsolv_util, is_iota_true) {
  auto vec = std::vector<int>{1, 2, 3, 4, 5};
  ASSERT_TRUE(is_iota(vec.begin(), vec.end(), 1));
}

TEST(itsolv_util, is_iota_false) {
  auto vec = std::vector<int>{1, 2, 3, 4, 5};
  ASSERT_FALSE(is_iota(vec.begin(), vec.end(), 0));
  auto vec2 = std::vector<int>{1, 3, 4, 5};
  ASSERT_FALSE(is_iota(vec2.begin(), vec2.end(), vec2[0]));
  auto vec3 = std::vector<int>{3, 2, 1};
  ASSERT_FALSE(is_iota(vec3.begin(), vec3.end(), vec3[0]));
}

TEST(itsolv_util, construct_zeroed_copy) {
  using R = std::vector<double>;
  using Q = std::list<double>;
  auto handler = molpro::linalg::array::create_default_handler<Q, R>();
  auto r = R(10);
  std::iota(begin(r), end(r), 0);
  auto q = construct_zeroed_copy(r, *handler);
  ASSERT_EQ(q.size(), r.size());
  ASSERT_THAT(q, ::testing::Each(::testing::DoubleEq(0)));
}

TEST(itsolv_util, delete_parameters__empty_params) {
  auto params = std::vector<int>{};
  ASSERT_NO_FATAL_FAILURE(delete_parameters({}, params));
  ASSERT_TRUE(params.empty());
}

TEST(itsolv_util, delete_parameters__empty_indices) {
  auto params = std::vector<int>{1, 2, 3, 4};
  auto params_ref = params;
  ASSERT_NO_FATAL_FAILURE(delete_parameters({}, params));
  ASSERT_THAT(params, ::testing::Pointwise(::testing::Eq(), params_ref));
}

TEST(itsolv_util, delete_parameters) {
  auto params = std::vector<int>{1, 2, 3, 4};
  auto indices = std::vector<int>{1, 3};
  auto params_ref = std::vector<int>{1, 3};
  ASSERT_NO_FATAL_FAILURE(delete_parameters(indices, params));
  ASSERT_THAT(params, ::testing::Pointwise(::testing::Eq(), params_ref));
}

TEST(StringFacet, toupper) {
  auto s = std::string("MixeD C@sE");
  auto sf = StringFacet{};
  auto s_upper = sf.toupper(s);
  ASSERT_EQ(s_upper, "MIXED C@SE");
}

TEST(StringFacet, tolower) {
  auto s = std::string("MixeD C@sE");
  auto sf = StringFacet{};
  auto s_lower = sf.tolower(s);
  ASSERT_EQ(s_lower, "mixed c@se");
}

TEST(StringFacet, crop_space) {
  auto core_text = std::string("some_words");
  for (std::string s : {" " + core_text, core_text + " ", " " + core_text + " "}) {
    StringFacet::crop_space(s);
    ASSERT_EQ(s, core_text);
  }
}