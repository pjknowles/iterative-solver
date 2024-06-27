#include <gtest/gtest.h>
#include <optional>

#include <molpro/linalg/itsolv/Interpolate.h>

using molpro::linalg::itsolv::Interpolate;

TEST(Interpolate, cubic_constructor) {
  // f = x(x-1/2)^2
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
  EXPECT_EQ(inter(0.5).f1, 0);
}

TEST(Interpolate, default_constructor) {
  // f = x(x-1/2)^2
  Interpolate::point p0{0, 0, 0.25};
  Interpolate::point p1{1, .25, 1.25};
  Interpolate inter(p0, p1);
  EXPECT_EQ(inter.interpolant(),"cubic");
  std::optional<std::string> empty;
  EXPECT_EQ(empty.has_value(),false);
  Interpolate inter_empty(p0,p1,empty);
  EXPECT_EQ(inter_empty.interpolant(),"cubic");
  std::optional<std::string> specified="cubic";
  EXPECT_EQ(specified.has_value(),true);
  Interpolate inter_specified(p0,p1,specified);
  EXPECT_EQ(inter_specified.interpolant(),"cubic");
}

TEST(Interpolate, minimize) {
  // f = x(x-1/2)^2
  for (const auto& interpolant : Interpolate::interpolants()) {
    if (interpolant != "morse") {
      Interpolate::point p0{0, 0, 0.25};
      Interpolate::point p1{1, .25, 1.25};
      Interpolate inter(p0, p1, interpolant, 0);
      EXPECT_NEAR(inter.minimize(0, 1, 5, 1000, false).x, 0.5, 1e-13);
      EXPECT_NEAR(inter.minimize(0, 1, 100, 100000, false).f, 0, 1e-13);
      EXPECT_NEAR(inter.minimize(0, 1, 100, 100000, false).f1, 0, 1e-13);
      EXPECT_GT(inter.minimize(0.3, 0.6, 100, 100000, false).f2, 0);
      EXPECT_NEAR(inter.minimize(0, -1, 100, 100000, false).x, -1, 1e-13);
      EXPECT_NEAR(inter.minimize(0.51, 1, 100, 100000, false).x, .51, 1e-13);
      EXPECT_NEAR(inter.minimize(0.1, 2, 100, 100000, false).x, .5, 1e-13);
      EXPECT_NEAR(inter.minimize(-1200, 2, 100, 100000, false).x, .5, 1e-13);
      EXPECT_NEAR(inter.minimize(.4, 200, 100, 100000, false).x, .5, 1e-13);
    }
  }
}

TEST(Interpolate, quadratic) {
  // f = (x-1/2)^2 + lambda x^3
  for (const auto& interpolant : Interpolate::interpolants()) {
    double lambda(1e-3);
//    std::cout << interpolant << std::endl;
    Interpolate::point p0{0, 0.25, -1};
//    Interpolate::point p0{-1, 2.25-lambda, -3+3*lambda};
    Interpolate::point p1{1, 0.25+lambda, 1+3*lambda};
    Interpolate inter(p0, p1, interpolant, 0);
//    std::cout << "found interpolant"<<inter<<std::endl;
    auto x0_expected = (-2.0+std::sqrt(4.0+12*lambda))/(6*lambda);
    auto minimum = inter.minimize(0, 1);
//    std::cout << "found minimum "<<minimum.x<<" value="<<minimum.f<< " slope="<<minimum.f1<<" second="<<minimum.f2<<std::endl;
    EXPECT_NEAR(minimum.x, x0_expected, 1e-8) << interpolant;
  }
}
