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

std::vector<std::string> Interpolate::interpolants() { return std::vector<std::string>{"cubic"}; }

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

Interpolate::point Interpolate::minimize(double xa, double xb, size_t bracket_grid, size_t max_bracket_grid) const {
  if (xa > xb)
    std::swap(xa, xb);
  //  std::cout << "Interpolate::minimize "<<xa<<" "<<xb<<", starting grid size "<<bracket_grid<<std::endl;
  for (size_t ngrid = bracket_grid; ngrid < std::max(bracket_grid, max_bracket_grid) + 1; ngrid *= 2) {
    auto gridstep = (xb - xa) / ngrid;
    auto plow = (*this)(xa);
    auto p0 = (*this)(xa).f > (*this)(xb).f ? plow : (*this)(xb);
    auto p1 = p0;
    for (int igrid = 0; igrid < ngrid; igrid++) {
      auto phigh = (*this)(plow.x + gridstep);
      //            std::cout << "Try bracket {"<<plow.x<<","<<plow.f<<","<<plow.f1<<"} :
      //            {"<<phigh.x<<","<<phigh.f<<","<<phigh.f1<<"}"<<std::endl;
      if (std::min(phigh.f, plow.f) < p0.f and plow.f1 <= 0 and phigh.f1 >= 0) {
        p1 = phigh;
        p0 = plow;
      }
      std::swap(plow, phigh);
    }
    //    std::cout << "Candidate bracket {"<<p0.x<<","<<p0.f<<","<<p0.f1<<"} :
    //    {"<<p1.x<<","<<p1.f<<","<<p1.f1<<"}"<<std::endl;
    if (p0.f1 < 0 and p1.f1 > 0) {
      //      std::cout << "Found bracket " << p0.x << " " << p1.x << std::endl;
      auto pnew = p1;
      auto tolerance = (std::nextafter(pnew.x, pnew.x + 1) - pnew.x) * 2;
      while (std::abs(p0.x - pnew.x) > tolerance) {
        //        std::cout << "iterate " << p0.x << ", " << p0.f << ", " << p0.f1 << ", " << p0.f2 << std::endl;
        pnew = (*this)((p1.x * p0.f1 - p0.x * p1.f1) / (p0.f1 - p1.f1));
        //        std::cout << "p0.x-p1.x " << p0.x - p1.x << ", pnew.x-p0.x" << pnew.x - p0.x << std::endl;
        //        std::cout << "p0 bound ranges " << std::nextafter(p0.x, p0.x + 1) - p0.x << ", "
        //                  << std::nextafter(p0.x, p0.x - 1) - p0.x << std::endl;
        if (pnew.f1 * p0.f1 < 0)
          std::swap(p0, p1);
        std::swap(p0, pnew);
      }
      //      std::cout << "return " << p0.x << ", " << p0.f << ", " << p0.f1 << ", " << p0.f2 << std::endl;
      return p0;
    }
  }
  // nothing found; return lowest end point
  return (*this)(xa).f > (*this)(xb).f ? (*this)(xb) : (*this)(xa);
}
