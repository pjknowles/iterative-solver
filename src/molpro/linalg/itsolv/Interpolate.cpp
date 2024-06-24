#include "Interpolate.h"
#include <cmath>
#include <iostream>
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/SolverFactory-implementation.h>
#include <molpro/linalg/itsolv/SolverFactory.h>

using namespace molpro::linalg::itsolv;

template <typename s>
std::ostream& operator<<(std::ostream& o, const std::vector<s>& v) {
  for (const auto& e : v)
    o << " " << e;
  return o;
}

using R = std::vector<double>;

Interpolate::point Morse(double y, const std::vector<double>& parameters) {
  Interpolate::point result;
  result.x = y;
  result.f = parameters[0] +
             (parameters[1] / 2) * std::pow((1 - std::exp(-parameters[2] * (y - parameters[3]))) / parameters[2], 2);
  result.f1 = (parameters[1] / parameters[2]) * std::exp(-parameters[2] * (y - parameters[3])) *
              (1 - std::exp(-parameters[2] * (y - parameters[3])));
  result.f2 = -parameters[1] * (1 - 2 * std::exp(-parameters[2] * (y - parameters[3])));
  return result;
}

class Morse_problem : public Problem<R> {
  Interpolate::point p0, p1;

public:
  Morse_problem(Interpolate::point p0, Interpolate::point p1) {
    this->p0 = p0;
    this->p1 = p1;
  }

  value_t residual(const R& parameters, R& residual) const override {
    auto pp0 = Morse(p0.x, parameters);
    auto pp1 = Morse(p1.x, parameters);
    residual[0] = pp0.f - p0.f;
    residual[1] = pp1.f - p1.f;
    residual[2] = pp0.f1 - p0.f1;
    residual[3] = pp1.f1 - p1.f1;
    //    std::cout << "residual, parameters " << parameters << std::endl;
    //    std::cout << "residual, residual " << residual << std::endl;
    return 0;
  }
};
Interpolate::Interpolate(point p0, point p1, std::string interpolant, int verbosity)
    : m_p0(std::move(p0)), m_p1(std::move(p1)), m_interpolant(std::move(interpolant)), m_parameters(4) {
  //  std::cout << "Interpolate constructor\n";
  //  std::cout << "Point " << p0.x << " " << p0.f << " " << p0.f1 << " " << p0.f2 << std::endl;
  //  std::cout << "Point " << p1.x << " " << p1.f << " " << p1.f1 << " " << p1.f2 << std::endl;
  if (m_interpolant == "cubic") {
    // c0 + c1(x-xbar) + c2(x-xbar)^2 + c3(x-xbar)^3 where xbar=(x0+x1)/2
    auto x1mx0 = m_p1.x - m_p0.x;
    auto f1pf0 = m_p1.f + m_p0.f;
    auto f1mf0 = m_p1.f - m_p0.f;
    auto g1pg0 = m_p1.f1 + m_p0.f1;
    auto g1mg0 = m_p1.f1 - m_p0.f1;
    m_parameters[0] = 0.5 * f1pf0 - 0.125 * g1mg0 * x1mx0;
    m_parameters[1] = -0.25 * g1pg0 + 1.5 * f1mf0 / x1mx0;
    m_parameters[2] = 0.5 * g1mg0 / x1mx0;
    m_parameters[3] = (-2 * f1mf0 + g1pg0 * x1mx0) / std::pow(x1mx0, 3);
  } else if (m_interpolant == "morse") {
    // L0 + (k/2a^2)*(1-exp(-a(y-y0)))^2
    // m_parameters: L0, k, a, y0
    R parameters(4);
    R residual(4);
    auto solver = molpro::linalg::itsolv::create_NonLinearEquations<R>("DIIS");
    auto cubic = Interpolate(p0, p1, "cubic", 0);
    auto cubic_minimum = cubic.minimize(p0.x, p1.x);
    auto cubic_at_minimum = cubic(cubic_minimum.x);
    //    parameters[0] = p0.f;
    //    parameters[1] = (p0.f1 - p1.f1) / (p0.x - p1.x);
    //    parameters[2] = 1e-10;
    //    parameters[3] = p0.x;
    //    std::cout << "cubic_at_minimum " << cubic_at_minimum.x << ", " << cubic_at_minimum.f << ", " <<
    //    cubic_at_minimum.f1
    //              << ", " << cubic_at_minimum.f2 << std::endl;
    //    std::cout << "cubic parameters " << cubic.m_parameters[0] << " " << cubic.m_parameters[1] << " "
    //              << cubic.m_parameters[2] << " " << cubic.m_parameters[3] << std::endl;
    m_parameters[1] = cubic_at_minimum.f2;
    m_parameters[2] = -3 * cubic.m_parameters[3] / (cubic_at_minimum.f2);
    m_parameters[3] = cubic_minimum.x;
    m_parameters[0] = cubic_at_minimum.f;
    //    std::fill(m_parameters.begin(), m_parameters.end(), double(1));
    Morse_problem problem(p0, p1);
    solver->set_verbosity(verbosity);
    if (!solver->solve(m_parameters, residual, problem))
      throw std::runtime_error("Cannot find Morse interpolant");
    solver->solution(m_parameters, residual);
  } else
    throw std::runtime_error("Unknown interpolant: " + m_interpolant);
}

std::vector<std::string> Interpolate::interpolants() { return std::vector<std::string>{"cubic", "morse"}; }

Interpolate::point Interpolate::operator()(double x) const {
  if (m_interpolant == "cubic") {
    auto xbar = 0.5 * (m_p1.x + m_p0.x);
    auto f = m_parameters[0] +
             (x - xbar) * (m_parameters[1] + (x - xbar) * (m_parameters[2] + (x - xbar) * m_parameters[3]));
    auto f1 = m_parameters[1] + (x - xbar) * (2 * m_parameters[2] + 3 * (x - xbar) * m_parameters[3]);
    auto f2 = 2 * m_parameters[2] + 6 * (x - xbar) * m_parameters[3];
    return point{x, f, f1, f2};
  } else if (m_interpolant == "morse") {
    return Morse(x, m_parameters);
  }
  throw std::logic_error("Unknown interpolant: " + m_interpolant);
}

Interpolate::point Interpolate::minimize_cubic() const {
//  std::cout << "minimize_cubic" << *this << std::endl;
  if (m_interpolant != "cubic")
    throw std::logic_error("minimize_cubic called with non-cubic interpolant");
  const auto c = m_parameters[1];
  const auto b = 2 * m_parameters[2];
  const auto a = 3 * m_parameters[3];
  auto discriminant = b*b / (4 * a*a) - c /  a;
//  std::cout << "a " << a << std::endl;
//  std::cout << "b " << b << std::endl;
//  std::cout << "c " << c << std::endl;
//  std::cout << "discriminant " << discriminant << std::endl;
  if (std::isnan(discriminant) || discriminant < 0)
    return {std::nan("unset")};
  auto xbar = 0.5 * (m_p1.x + m_p0.x);
  Interpolate::point pm = (*this)(xbar - (b / (2 * a)) + std::sqrt(discriminant));
  Interpolate::point pp = (*this)(xbar - (b / (2 * a)) - std::sqrt(discriminant));
//  std::cout << "extrema \n" << pm << "\n" << pp << std::endl;
  return pm.f < pp.f ? pm : pp;
}

Interpolate::point Interpolate::minimize(double xa, double xb, size_t bracket_grid, size_t max_bracket_grid, bool analytic) const {
  if (xa > xb)
    std::swap(xa, xb);
//  std::cout << "Interpolate::minimize " << xa << " " << xb << ", starting grid size " << bracket_grid << std::endl;
//  std::cout << "Interpolant:" << *this << std::endl;
  if (analytic && m_interpolant == "cubic")
    return minimize_cubic();
//    std::cout << "minimize_cubic() " << minimize_cubic() << std::endl;
  for (size_t ngrid = bracket_grid; ngrid < std::max(bracket_grid, max_bracket_grid) + 1; ngrid *= 2) {
    auto gridstep = (xb - xa) / ngrid;
    auto plow = (*this)(xa);
    //    std::cout << "plow " << plow << std::endl;
    auto p0 = (*this)(xa).f > (*this)(xb).f ? plow : (*this)(xb);
    auto p1 = p0;
    for (size_t igrid = 0; igrid < ngrid; igrid++) {
      auto phigh = (*this)(plow.x + gridstep);
      //      std::cout << "Try bracket {" << plow.x << "," << plow.f << "," << plow.f1 << "} : {" << phigh.x << "," <<
      //      phigh.f
      //                << "," << phigh.f1 << "}" << std::endl;
      if (std::min(phigh.f, plow.f) < p0.f and plow.f1 <= 0 and phigh.f1 >= 0) {
        p1 = phigh;
        p0 = plow;
      }
      std::swap(plow, phigh);
    }
    //    std::cout << "Candidate bracket {" << p0.x << "," << p0.f << "," << p0.f1 << "} : {" << p1.x << "," << p1.f <<
    //    ","
    //              << p1.f1 << "}" << std::endl;
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

std::ostream& molpro::linalg::itsolv::operator<<(std::ostream& os, const Interpolate& interpolant) {
  for (const auto& parameter : interpolant.parameters())
    os << " " << parameter;
  return os;
}
std::ostream& molpro::linalg::itsolv::operator<<(std::ostream& os, const Interpolate::point& p) {
  os << "x=" << p.x << ", value=" << p.f << ", gradient=" << p.f1 << ", curvature=" << p.f2;
  return os;
}