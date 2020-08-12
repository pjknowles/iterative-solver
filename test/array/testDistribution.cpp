#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/array/util/Distribution.h>

using ::testing::DoubleEq;
using ::testing::Pointwise;

using molpro::linalg::array::util::Distribution;
using molpro::linalg::array::util::make_distribution_spread_remainder;

TEST(Distribution, constructor) {
  auto chunk_borders = std::vector<size_t>{0, 3, 7, 9};
  auto d = Distribution<size_t>(chunk_borders);
  ASSERT_THAT(d.chunk_borders(), Pointwise(DoubleEq(), chunk_borders));
}

TEST(Distribution, size) {
  auto chunk_borders = std::vector<size_t>{0, 3, 7, 9};
  auto d = Distribution<size_t>(chunk_borders);
  ASSERT_EQ(d.size(), chunk_borders.size() - 1);
}

TEST(Distribution, border) {
  auto chunk_borders = std::vector<size_t>{0, 3, 7, 9};
  auto d = Distribution<size_t>(chunk_borders);
  ASSERT_TRUE((d.border() == std::pair<size_t, size_t>{0, 9}));
}

TEST(Distribution, range) {
  auto chunk_borders = std::vector<size_t>{0, 3, 7, 9};
  auto d = Distribution<size_t>(chunk_borders);
  EXPECT_EQ(d.range(0), (std::pair<size_t, size_t>{0, 3}));
  EXPECT_EQ(d.range(1), (std::pair<size_t, size_t>{3, 7}));
  EXPECT_EQ(d.range(2), (std::pair<size_t, size_t>{7, 9}));
}

TEST(Distribution, cover) {
  auto chunk_borders = std::vector<size_t>{0, 3, 7, 9};
  auto d = Distribution<size_t>(chunk_borders);
  auto ref_cover = std::vector<int>{0, 0, 0, 1, 1, 1, 1, 2, 2};
  auto cover = std::vector<int>(ref_cover.size());
  for (size_t i = 0; i < cover.size(); ++i)
    cover[i] = d.cover(i);
  ASSERT_THAT(cover, Pointwise(DoubleEq(), ref_cover));
  ASSERT_EQ(d.cover(10), d.size());
}

TEST(Distribution, cover_2params) {
  auto chunk_borders = std::vector<size_t>{0, 3, 7, 9};
  auto d = Distribution<size_t>(chunk_borders);
  auto ref_cover = std::vector<int>{0, 0, 0, 1, 1, 1, 1, 2, 2};
  auto cover = std::vector<int>(ref_cover.size());
  for (size_t i = 0; i < cover.size(); ++i)
    for (size_t j = i; j < cover.size(); ++j)
      ASSERT_TRUE((d.cover(i, j) == std::pair<int, int>{ref_cover[i], ref_cover[j]}));
}

TEST(Distribution, make_distribution_spread_remainder) {
  auto d = make_distribution_spread_remainder<size_t>(11, 3);
  auto ref_chunk_borders = std::vector<size_t>{0, 4, 8, 11};
  ASSERT_THAT(d.chunk_borders(), Pointwise(DoubleEq(), ref_chunk_borders));
}
