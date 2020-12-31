#include <gtest/gtest.h>

#include <molpro/linalg/itsolv/Interpolate.h>

using molpro::linalg::itsolv::Interpolate;

TEST(Interpolate, cubic_constructor) {
  Interpolate::point p0{0, 0, 0.25};
  Interpolate::point p1{1, .25, 1.25};
  Interpolate inter(p0, p1, "cubic");
  EXPECT_EQ(inter(p0.x).x, p0.x);
  EXPECT_EQ(inter(p0.x).f, p0.f);
  EXPECT_EQ(inter(p0.x).f1, p0.f1);
  EXPECT_EQ(inter(p0.x).f2, -2);
  EXPECT_EQ(inter(p1.x).x, p1.x);
  EXPECT_EQ(inter(p1.x).f, p1.f);
  EXPECT_EQ(inter(p1.x).f1, p1.f1);
  EXPECT_EQ(inter(p1.x).f2, 4);
  EXPECT_EQ(inter(0.5).f1,0);
}

TEST(Interpolate, minimize) {
  Interpolate::point p0{0, 0, 0.25};
  Interpolate::point p1{1, .25, 1.25};
  Interpolate inter(p0, p1, "cubic");
  EXPECT_NEAR(inter.minimize(0,1).x,0.5,1e-8);
  EXPECT_NEAR(inter.minimize(0,1).f1,0,1e-8);
  EXPECT_GT(inter.minimize(0,1).f2,0);
}
