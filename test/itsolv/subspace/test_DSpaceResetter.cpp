#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/itsolv/DSpaceResetter.h>

using molpro::linalg::itsolv::detail::DoReset;

TEST(itsolv_reset_dspace, DoReset_constructor) {
  auto d = DoReset();
  ASSERT_FALSE(d.value());
  const size_t n = 7;
  auto d2 = DoReset(n);
  ASSERT_FALSE(d2.value());
}

TEST(itsolv_reset_dspace, DoReset_set_get_nreset) {
  const size_t n = 7, m = 3;
  auto d = DoReset(n);
  ASSERT_EQ(d.get_nreset(), n);
  d.set_nreset(m);
  ASSERT_EQ(d.get_nreset(), m);
}

TEST(itsolv_reset_dspace, DoReset_update) {
  const size_t n = 2;
  const unsigned int init_max_q = 10, nD = 3;
  unsigned int max_q = init_max_q;
  auto d = DoReset(n);
  d.update(0, max_q, nD);
  ASSERT_FALSE(d.value());
  ASSERT_EQ(max_q, init_max_q);
  d.update(1, max_q, nD);
  ASSERT_TRUE(d.value());
  ASSERT_EQ(max_q, init_max_q + nD);
  d.update(2, max_q, nD);
  ASSERT_TRUE(d.value());
  ASSERT_EQ(max_q, init_max_q + nD);
  d.update(3, max_q, 0);
  ASSERT_FALSE(d.value());
  ASSERT_EQ(max_q, init_max_q);
}
