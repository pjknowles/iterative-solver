#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/itsolv/options_map.h>
#include <molpro/linalg/itsolv/util.h>

using molpro::linalg::itsolv::options_map;
using molpro::linalg::itsolv::util::capitalize_keys;
using molpro::linalg::itsolv::util::StringFacet;

TEST(options_map, capitalize_keys) {
  auto sf = StringFacet{};
  auto opt = options_map{{"key1", "val1"}, {"kEy2", "Val2"}, {"KEY3", "VAL3"}};
  auto opt_ref = options_map{{"KEY1", "val1"}, {"KEY2", "Val2"}, {"KEY3", "VAL3"}};
  auto opt_upper = capitalize_keys(opt, sf);
  ASSERT_EQ(opt_upper, opt_ref);
}
