#include "Interpolate.h"
#include <cmath>
#include <iostream>

using namespace molpro::linalg::itsolv;
Interpolate::Interpolate(point p0, point p1, std::string interpolant)
    : m_p0(std::move(p0)), m_p1(std::move(p1)), m_interpolant(std::move(interpolant)) {
  if (m_interpolant == "cubic") {
    // c0 + c1(x-xbar) + c2(x-xbar)^2 + c3(x-xbar)^3 where xbar=(x0+x1)/2
    auto x1mx0 = m_p1.x - m_p0.x;
    auto f1pf0 = m_p1.f + m_p0.f;
    auto f1mf0 = m_p1.f - m_p0.f;
    auto g1pg0 = m_p1.f1 + m_p0.f1;
    auto g1mg0 = m_p1.f1 - m_p0.f1;
    m_c0 = 0.5 * f1pf0 - 0.125 * g1mg0 * x1mx0;
    m_c1 = -0.25 * g1pg0 + 1.5 * f1mf0 / x1mx0;
    m_c2 = 0.5 * g1mg0 / x1mx0;
    m_c3 = (-2 * f1mf0 + g1pg0 * x1mx0) / std::pow(x1mx0, 3);
  } else
    throw std::runtime_error("Unknown interpolant: " + m_interpolant);
}

Interpolate::point Interpolate::operator()(double x) const {
  if (m_interpolant == "cubic") {
    auto xbar = 0.5 * (m_p1.x + m_p0.x);
    auto f = m_c0 + (x - xbar) * (m_c1 + (x - xbar) * (m_c2 + (x - xbar) * m_c3));
    auto f1 = m_c1 + (x - xbar) * (2 * m_c2 + 3 * (x - xbar) * m_c3);
    auto f2 = 2 * m_c2 + 6 * (x - xbar) * m_c3;
    return point{x, f, f1, f2};
  }
  throw std::logic_error("Unknown interpolant: " + m_interpolant);
}

Interpolate::point Interpolate::minimize(double xmin, double xmax) const {
  double x = xmin;
  return (*this)(x);
}
