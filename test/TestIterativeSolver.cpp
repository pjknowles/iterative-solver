#include "test.h"
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "molpro/linalg/IterativeSolver.h"
#include "molpro/linalg/PagedArray.h"
#include "molpro/linalg/SimpleArray.h"
#include "molpro/linalg/array/ArrayHandlerIterable.h"
#include <Eigen/Dense>
#include <cmath>
#include <limits>
#include <molpro/iostream.h>
#include <numeric>
#include <regex>
#include <vector>

TEST(TestIterativeSolver, small_eigenproblem) {
  for (size_t n = 1; n < 9; n++) {
    for (size_t nroot = 1; nroot <= n && nroot < 3; nroot++) {
      Eigen::MatrixXd m(n, n);
      for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++)
          m(i, j) = 1 + (i + j) * std::sqrt(double(i + j));
        m(i, i) += i;
      }
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> denseSolver(m);
      auto val = denseSolver.eigenvalues();

      molpro::linalg::SimpleArray<double> mm(n);
      std::vector<molpro::linalg::SimpleArray<double>> x, g;
      molpro::linalg::array::ArrayHandlerIterable<molpro::linalg::SimpleArray<double>> handler{};
      molpro::linalg::LinearEigensystem<molpro::linalg::SimpleArray<double>> solver;
      solver.m_verbosity = 1;
      solver.setThresholds(1e-12);
      if (solver.m_verbosity > 0)
        std::cout << "Test n=" << n << ", nroot=" << nroot << std::endl;
      if (solver.m_verbosity > 1)
        std::cout << "Matrix\n" << m << std::endl;
      for (size_t root = 0; root < nroot; root++) {
        x.emplace_back(n);
        x.back().scal(0);
        x.back()[root] = 1;
        g.emplace_back(n);
      }
      size_t nwork = nroot;
      for (size_t iter = 0; iter < n + 10; iter++) {
        if (solver.m_verbosity > 1)
          std::cout << "start iteration " << iter << ", nwork=" << nwork << std::endl;
        for (size_t root = 0; root < nwork; root++) {
          g[root].scal(0);
          for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
              g[root][j] += m(j, i) * x[root][i];
        }

        //        std::cout << "eigenvector "<<0<<active[0]<<" before addVector"; for (size_t i = 0; i < n; i++)
        //        std::cout << " " << x[0][i]; std::cout << std::endl;
        auto nwork = solver.addVector(x, g);
        if (nwork == 0)
          break;
        for (size_t root = 0; root < nwork; root++) {
          if (solver.m_verbosity > 1) {
            std::cout << "eigenvector " << root << " before update";
            for (size_t i = 0; i < n; i++)
              std::cout << " " << x[root][i];
            std::cout << std::endl;
          }
          {
            for (size_t i = 0; i < n; i++)
              x[root][i] -= g[root][i] / (m(i, i) - solver.eigenvalues()[root] + 1e-13);
            if (solver.m_verbosity > 2) {
              std::cout << "residual " << root << " ";
              for (size_t i = 0; i < n; i++)
                std::cout << " " << g[root][i];
              std::cout << std::endl;
            }
            if (solver.m_verbosity > 1) {
              std::cout << "eigenvector " << root << " ";
              for (size_t i = 0; i < n; i++)
                std::cout << " " << x[root][i];
              std::cout << std::endl;
            }
          }
        }
        //                std::cout << "eigenvector "<<0<<" before endIteration"; for (size_t i = 0; i < n; i++)
        //                std::cout << " " << x[0][i]; std::cout << std::endl;
        //                auto conv = (solver.endIteration(x, g, active));
        //        std::cout << "eigenvector "<<0<<active[0]<<" after endIteration"; for (size_t i = 0; i < n; i++)
        //        std::cout << " " << x[0][i]; std::cout << std::endl; if (conv) break;
        //        if (*std::max_element(solver.errors().begin(),solver.errors().end()) < solver.m_thresh and
        //        solver.endIteration(x, g))
      }
      //  std::cout << "Error={ "; for (const auto& e : solver.errors()) std::cout << e << " "; std::cout << "} after "
      //  << solver.iterations() << " iterations" << std::endl; std::cout << "Actual eigenvalues\n"<<val<<std::endl;
      //      EXPECT_THAT(active,
      //                  ::testing::Pointwise(::testing::Eq(), std::vector<bool>(nroot, false)));
      EXPECT_THAT(solver.errors(),
                  ::testing::Pointwise(::testing::DoubleNear(1e-10), std::vector<double>(nroot, double(0))));
      EXPECT_THAT(solver.eigenvalues(), ::testing::Pointwise(::testing::DoubleNear(1e-10),
                                                             std::vector<double>(val.data(), val.data() + nroot)));
      std::vector<int> roots(solver.m_roots);
      std::iota(roots.begin(), roots.end(), 0);
      solver.solution(roots, x, g);
      for (size_t root = 0; root < solver.m_roots; root++) {
        if (solver.m_verbosity > 2) {
          std::cout << "eigenvector " << root << " eigenvalue=" << solver.eigenvalues()[root]
                    << " converged=" << solver.errors()[root] << ":";
          for (size_t i = 0; i < n; i++)
            std::cout << " " << x[root][i];
          std::cout << std::endl;
        }
        std::vector<double> r(n);
        g[root].get(r.data(), n, 0);
        EXPECT_THAT(r, ::testing::Pointwise(::testing::DoubleNear(1e-5), std::vector<double>(n, double(0))));
        if (solver.m_verbosity > 1)
          for (size_t soot = 0; soot <= root; soot++)
            std::cout << "Eigenvector overlap " << root << " " << soot << " " << handler.dot(x[root],x[soot]) << std::endl;
        for (size_t soot = 0; soot < root; soot++)
          EXPECT_LE(std::abs(handler.dot(x[root],x[soot])),
                    1e-8); // can't expect exact orthogonality when last thing might have been an update
        EXPECT_THAT(std::abs(handler.dot(x[root],x[root])), ::testing::DoubleNear(1, 1e-10));
      }
    }
  }
}

TEST(TestIterativeSolver, small_nonhermitian_eigenproblem) {
  for (size_t n = 2; n < 8; n++) {
    for (size_t nroot = 1; nroot <= n && nroot < 6; nroot++) {
      Eigen::MatrixXd m(n, n);
      for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++)
          m(i, j) = 2 * i + j + 3;
      }
      Eigen::EigenSolver<Eigen::MatrixXd> denseSolver(m);
      auto val = denseSolver.eigenvalues();

      molpro::linalg::SimpleArray<double> mm(n);
      std::vector<molpro::linalg::SimpleArray<double>> x, g;
      molpro::linalg::LinearEigensystem<molpro::linalg::SimpleArray<double>> solver;
      solver.m_verbosity = 1;
      solver.setThresholds(1e-13);
      if (solver.m_verbosity > 0)
        std::cout << "Test n=" << n << ", nroot=" << nroot << std::endl;
      for (size_t root = 0; root < nroot; root++) {
        x.emplace_back(n);
        x.back().scal(0);
        x.back()[root] = 1;
        g.emplace_back(n);
      }
      for (size_t iter = 0; iter < n + 10; iter++) {
        for (size_t root = 0; root < x.size(); root++) {
          g[root].scal(0);
          for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
              g[root][j] += m(j, i) * x[root][i];
        }

        //        std::cout << "eigenvector "<<0<<" before addVector"; for (size_t i = 0; i < n; i++) std::cout << " "
        //        << x[0][i]; std::cout << std::endl; std::cout << "action "<<0<<" before addVector"; for (size_t i = 0;
        //        i < n; i++) std::cout << " " << g[0][i]; std::cout << std::endl; if (x.size() > 1) {
        //          std::cout << "eigenvector "<<1<<" before addVector"; for (size_t i = 0; i < n; i++) std::cout << " "
        //          << x[1][i]; std::cout << std::endl;
        //        std::cout << "action "<<1<<" before addVector"; for (size_t i = 0; i < n; i++) std::cout << " " <<
        //        g[1][i]; std::cout << std::endl;
        //        }
        auto nwork = solver.addVector(x, g);
        for (size_t root = 0; root < nwork; root++) {
          if (solver.m_verbosity > 1) {
            std::cout << "eigenvector " << root << " before update";
            for (size_t i = 0; i < n; i++)
              std::cout << " " << x[root][i];
            std::cout << std::endl;
          }
          for (size_t i = 0; i < n; i++)
            x[root][i] -= g[root][i] / (m(i, i) - solver.eigenvalues()[root] + 1e-13);
          if (solver.m_verbosity > 2) {
            std::cout << "residual " << root << " ";
            for (size_t i = 0; i < n; i++)
              std::cout << " " << g[root][i];
            std::cout << std::endl;
          }
          if (solver.m_verbosity > 1) {
            std::cout << "eigenvector " << root << " ";
            for (size_t i = 0; i < n; i++)
              std::cout << " " << x[root][i];
            std::cout << std::endl;
          }
        }
        //        std::cout << "eigenvector "<<0<<active[0]<<" before endIteration"; for (size_t i = 0; i < n; i++)
        //        std::cout << " " << x[0][i]; std::cout << std::endl; auto conv = (solver.endIteration(x, g, active));
        //        std::cout << "eigenvector "<<0<<active[0]<<" after endIteration"; for (size_t i = 0; i < n; i++)
        //        std::cout << " " << x[0][i]; std::cout << std::endl; if (conv) break;
        if (nwork == 0)
          break;
      }
      //  std::cout << "Error={ "; for (const auto& e : solver.errors()) std::cout << e << " "; std::cout << "} after "
      //  << solver.iterations() << " iterations" << std::endl; std::cout << "Actual eigenvalues\n"<<val<<std::endl;
      //      EXPECT_THAT(active,
      //                  ::testing::Pointwise(::testing::Eq(), std::vector<bool>(nroot, false)));
      EXPECT_THAT(solver.errors(),
                  ::testing::Pointwise(::testing::DoubleNear(1e-10), std::vector<double>(nroot, double(0))));
      //      EXPECT_THAT(solver.eigenvalues(),
      //                  ::testing::Pointwise(::testing::DoubleNear(1e-10),
      //                                       std::vector<double>(val.real().data(), val.real().data() + nroot)));
      std::vector<int> roots;
      for (auto root = 0; root < solver.m_roots; root++)
        roots.push_back(root);
      solver.solution(roots, x, g);
      for (size_t root = 0; root < solver.m_roots; root++) {
        if (solver.m_verbosity > 1) {
          std::cout << "eigenvector " << root << " converged=" << solver.errors()[root] << ":";
          for (size_t i = 0; i < n; i++)
            std::cout << " " << x[root][i];
          std::cout << std::endl;
          std::cout << "residual:";
          for (size_t i = 0; i < n; i++)
            std::cout << " " << g[root][i];
          std::cout << std::endl;
        }
      }
    }
  }
}

TEST(TestIterativeSolver, linear_equations) {
  for (size_t n = 1; n < 9; n++) {
    for (size_t nroot = 1; nroot <= n && nroot < 10; nroot++) {
      Eigen::MatrixXd m(n, n);
      for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
          m(i, j) = 1 + (i + j) * std::sqrt(double(i + j));
      for (size_t i = 0; i < n; i++)
        m(i, i) += i;

      molpro::linalg::SimpleArray<double> mm(n);
      std::vector<molpro::linalg::SimpleArray<double>> x, g, rhs;
      for (size_t root = 0; root < nroot; root++) {
        x.emplace_back(n);
        x.back().scal(0);
        x.back()[root] = 1;
        g.emplace_back(n);
        rhs.emplace_back(n);
        rhs.back().scal(0);
        rhs.back()[root] = 1;
        Eigen::VectorXd erhs(n);
        rhs[root].get(&erhs(0), n, 0);
        auto trueSolution = m.colPivHouseholderQr().solve(erhs).eval();
        rhs.back()[root] = 1 / trueSolution(root);
      }
      molpro::linalg::LinearEquations<molpro::linalg::SimpleArray<double>> solver(rhs);
      solver.m_verbosity = 1;
      solver.setThresholds(1e-13);
      if (solver.m_verbosity > 0)
        std::cout << "Test n=" << n << ", nroot=" << nroot << std::endl;
      if (solver.m_verbosity > 1)
        std::cout << "Matrix:\n" << m << std::endl;
      for (size_t iter = 0; iter < n + 5; iter++) {
        for (size_t root = 0; root < x.size(); root++) {
          g[root].scal(0);
          for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
              g[root][j] += m(j, i) * x[root][i];
        }

        //        std::cout << "solution "<<0<<solver.active()[0]<<" before addVector"; for (size_t i = 0; i < n; i++)
        //        std::cout << " " << x[0][i]; std::cout << std::endl;
        auto nwork = solver.addVector(x, g);
        for (const auto& root : solver.working_set()) {
          if (solver.m_verbosity > 1) {
            std::cout << "solution " << root << " before update";
            for (size_t i = 0; i < n; i++)
              std::cout << " " << x[root][i];
            std::cout << std::endl;
          }
          for (size_t i = 0; i < n; i++)
            x[root][i] -= g[root][i] / (m(i, i) //                  - solver.eigenvalues()[root]
                                        + 1e-13);
          if (solver.m_verbosity > 2) {
            std::cout << "residual " << root << " ";
            for (size_t i = 0; i < n; i++)
              std::cout << " " << g[root][i];
            std::cout << std::endl;
          }
          if (solver.m_verbosity > 1) {
            std::cout << "solution " << root << " ";
            for (size_t i = 0; i < n; i++)
              std::cout << " " << x[root][i];
            std::cout << std::endl;
          }
        }
        //        std::cout << "eigenvector "<<0<<active[0]<<" before endIteration"; for (size_t i = 0; i < n; i++)
        //        std::cout << " " << x[0][i]; std::cout << std::endl; auto conv = (solver.endIteration(x, g, active));
        //        std::cout << "eigenvector "<<0<<active[0]<<" after endIteration"; for (size_t i = 0; i < n; i++)
        //        std::cout << " " << x[0][i]; std::cout << std::endl; if (conv) break;
        if (nwork == 0)
          break;
      }
      //      EXPECT_THAT(solver.active(),
      //                  ::testing::Pointwise(::testing::Eq(), std::vector<bool>(nroot, false)));
      EXPECT_THAT(solver.errors(),
                  ::testing::Pointwise(::testing::DoubleNear(1e-10), std::vector<double>(nroot, double(0))));
      //                                       std::vector<double>(val.data(), val.data() + nroot)));
      std::vector<int> roots;
      for (size_t root = 0; root < solver.m_roots; root++)
        roots.push_back(root);
      solver.solution(roots, x, g);
      for (size_t root = 0; root < solver.m_roots; root++) {
        if (solver.m_verbosity > 1) {
          std::cout << "solution " << root << " converged=" << solver.errors()[root] << ":";
          for (size_t i = 0; i < n; i++)
            std::cout << " " << x[root][i];
          std::cout << std::endl;
        }
        Eigen::VectorXd erhs(n);
        rhs[root].get(&erhs(0), n, 0);
        auto trueSolution = m.colPivHouseholderQr().solve(erhs).eval();
        //        std::cout << "RHS\n"<<erhs<<std::endl;
        //        std::cout << "trueSolution\n"<<trueSolution<<std::endl;
        std::vector<double> r(n);
        x[root].get(r.data(), n, 0);
        EXPECT_THAT(r, ::testing::Pointwise(::testing::DoubleNear(1e-5),
                                            std::vector<double>(&trueSolution(0), &trueSolution(0) + n)));
      }
    }
  }
}

#include <regex>
class RosenbrockTest {
public:
  using ptype = molpro::linalg::SimpleArray<double>;
  using scalar = typename molpro::linalg::Optimize<ptype>::scalar_type;
  static constexpr double Rosenbrock_a = 1;
  static constexpr double Rosenbrock_b = 1;
  std::string algorithm;
  struct {
    scalar operator()(const ptype& psx, ptype& outputs) const {
      size_t n = 2;
      std::vector<scalar> psxk(n);
      std::vector<scalar> output(n);

      psx.get(&(psxk[0]), n, 0);
      output[0] = (2 * psxk[0] - 2 * Rosenbrock_a + 4 * Rosenbrock_b * psxk[0] * (psxk[0] * psxk[0] - psxk[1]));
      output[1] = (2 * Rosenbrock_b * (psxk[1] - psxk[0] * psxk[0])); // Rosenbrock
      outputs.put(&(output[0]), n, 0);
      return (Rosenbrock_a - psxk[0]) * (Rosenbrock_a - psxk[0]) +
             Rosenbrock_b * (psxk[1] - psxk[0] * psxk[0]) * (psxk[1] - psxk[0] * psxk[0]);
    }
  } _Rosenbrock_residual;

  struct {
    bool approximateHessian = true;
    void operator()(ptype& psc, const ptype& psg, std::vector<scalar> shift, bool append = true) const {
      size_t n = 2;
      std::vector<scalar> psck(n);
      std::vector<scalar> psgk(n);
      psg.get(&psgk[0], n, 0);
      if (append)
        psc.get(&psck[0], n, 0);
      else
        psck.assign(0, 0);
      //      molpro::cout << "Rosenbrock updater, initial psc=" << psc << std::endl;
      //      molpro::cout << "Rosenbrock updater, initial psg=" << psg << std::endl;
      //      molpro::cout << "Rosenbrock updater, approximateHessian=" << approximateHessian << std::endl;
      const auto& x = psck[0];
      const auto& y = psck[1];
      scalar dx, dy;
      if (approximateHessian) {
        //      dx = -psgk[0] / 7*Rosenbrock_b;
        //      dy = -psgk[1] / 2*Rosenbrock_b;
        //        dx =
        //            -psgk[0] / (4 + 8 * Rosenbrock_b * (x * x - y));
        //        dy =
        //            -psgk[1] / (4 * Rosenbrock_b)
        //                - psgk[1] * x * x / (1 + 2 * Rosenbrock_b * (x * x - y));
        dx = dy = (-psgk[0] - psgk[1]) / 4;
      } else {
        dx = -psgk[0] / (4 + 8 * Rosenbrock_b * (x * x - y)) - psgk[1] / (4 + 4 * Rosenbrock_b * (x * x - y));
        dy = -psgk[1] / (4 * Rosenbrock_b) - psgk[1] * x * x / (1 + 2 * Rosenbrock_b * (x * x - y)) -
             psgk[0] / (4 + 4 * Rosenbrock_b * (x * x - y));
      }
      psck[0] += dx;
      psck[1] += dy;
      psc.put(&psck[0], n, 0);
      //      molpro::cout << "Rosenbrock updater, new psc=" << psc << std::endl;
    }
  } _Rosenbrock_updater;
  bool run(std::string method) {

    ptype x(2);
    ptype g(2);
    ptype hg(2);
    const double difficulty = .9;
    const int verbosity = 1;

    if (verbosity >= 0)
      molpro::cout << "Test Optimize, method=" << method << ", difficulty=" << difficulty << std::endl;
    molpro::linalg::Optimize<ptype> d(regex_replace(method, std::regex("-.*"), ""));
    d.m_verbosity = verbosity - 1;
    d.m_options["convergence"] = "residual";
    std::vector<scalar> xxx(2);
    xxx[0] = xxx[1] = 1 - difficulty; // initial guess
    x.put(&xxx[0], 2, 0);
    //    molpro::cout << "initial guess " << x << std::endl;
    _Rosenbrock_updater.approximateHessian = true;
    bool converged = false;
    for (int iteration = 1; iteration < 50 && not converged; iteration++) {
      auto value = _Rosenbrock_residual(x, g);
      std::vector<scalar> shift;
      shift.push_back(1e-10);
      if (d.addValue(x, value, g))
        _Rosenbrock_updater(x, g, shift);
      converged = d.endIteration(x, g);
      x.get(&xxx[0], 2, 0);
      if (verbosity >= 0)
        molpro::cout << "iteration " << iteration << ", Distance from solution = "
                     << std::sqrt((xxx[0] - Rosenbrock_a) * (xxx[0] - Rosenbrock_a) +
                                  (xxx[1] - Rosenbrock_a) * (xxx[1] - Rosenbrock_a))
                     << ", error = " << d.errors().front() << ", converged? " << converged << ", value= " << value
                     << std::endl;
    }

    x.get(&xxx[0], 2, 0);
    molpro::cout << "Distance from solution = "
                 << std::sqrt((xxx[0] - Rosenbrock_a) * (xxx[0] - Rosenbrock_a) +
                              (xxx[1] - Rosenbrock_a) * (xxx[1] - Rosenbrock_a))
                 << std::endl;

    return converged;
  }
};

TEST(Rosenbrock_BFGS, DISABLED_Optimize) { ASSERT_TRUE(RosenbrockTest().run("L-BFGS")); }
TEST(Rosenbrock_null, DISABLED_Optimize) { ASSERT_TRUE(RosenbrockTest().run("null")); }

template <class T>
std::ostream& operator<<(std::ostream& o, const molpro::linalg::SimpleArray<T>& a) {
  std::vector<T> v;
  v.resize(a.size());
  a.get(v.data(), v.size(), 0);
  o << v;
  return o;
}

class MonomialTest {
public:
  using ptype = molpro::linalg::SimpleArray<double>;
  using scalar = typename molpro::linalg::Optimize<ptype>::scalar_type;
  struct {
    double power;
    double normPower;
    scalar operator()(const ptype& psx, ptype& outputs) const {
      size_t n = psx.size();
      std::vector<scalar> psxk(n);
      std::vector<scalar> output(n);

      psx.get(&(psxk[0]), n, 0);
      double fp = 0;
      for (size_t i = 0; i < n; i++)
        fp += (i + 1) * std::pow(psxk[i], power);
      for (size_t i = 0; i < n; i++) {
        output[i] = psxk[i] == 0 ? 0
                                 : power * (i + 1) * std::pow(psxk[i], power - 1) * (normPower / power) *
                                       pow(fp, normPower / power - 1);
      }
      outputs.put(&(output[0]), n, 0);
      return fp;
    }
  } _Monomial_residual;

  struct {
    void operator()(ptype& psc, const ptype& psg, std::vector<scalar> shift, bool append = true) const {
      size_t n = psc.size();
      std::vector<scalar> psck(n);
      std::vector<scalar> psgk(n);
      psg.get(&psgk[0], n, 0);
      if (append)
        psc.get(&psck[0], n, 0);
      else
        psck.assign(0, 0);
      for (size_t i = 0; i < n; i++)
        psck[i] -= psgk[i];
      psc.put(&psck[0], n, 0);
    }
  } _Monomial_updater;
  bool run(double power, double normPower) {
    _Monomial_residual.power = power;
    _Monomial_residual.normPower = normPower;
    constexpr size_t n = 1;
    ptype x(n);
    ptype g(n);
    ptype hg(n);
    const int verbosity = 1;

    molpro::linalg::Optimize<ptype> d("L-BFGS");
    //    IterativeSolver::DIIS<ptype> d;
    d.m_verbosity = verbosity - 1;
    d.m_options["convergence"] = "residual";
    d.m_Wolfe_1 = 1e20;
    d.m_linesearch_tolerance = 0;
    std::vector<scalar> xxx(n, .1);
    x.put(&xxx[0], n, 0);
    if (verbosity > 1)
      molpro::cout << "initial guess " << x << std::endl;
    bool converged = false;
    for (int iteration = 1; iteration < 100 && not converged; iteration++) {
      auto value = _Monomial_residual(x, g);
      std::vector<scalar> shift;
      shift.push_back(1e-10);
      if (verbosity > 1)
        molpro::cout << "x" << x << std::endl;
      if (verbosity > 1)
        molpro::cout << "g" << g << std::endl;
      if (d.addValue(x, value, g))
        _Monomial_updater(x, g, shift);
      converged = d.endIteration(x, g);
      x.get(&xxx[0], n, 0);
      if (verbosity >= 0)
        molpro::cout << "iteration " << iteration << ", error = " << d.errors().front() << ", converged? " << converged
                     << ", value= " << value << std::endl;
    }

    //    x.get(&xxx[0], n, 0);
    //    molpro::cout << "final solution:" ;
    //    for (const auto& xxxx : xxx) molpro::cout << " "<<xxxx; molpro::cout <<std::endl;

    return converged;
  }
};

TEST(Monomial_22, Optimize) { ASSERT_TRUE(MonomialTest().run(2, 2)); }
TEST(Monomial_44, Optimize) { ASSERT_TRUE(MonomialTest().run(4, 4)); }
TEST(Monomial_42, Optimize) { ASSERT_TRUE(MonomialTest().run(4, 2)); }

class MorseTest {
public:
  using ptype = molpro::linalg::SimpleArray<double>;
  using scalar = typename molpro::linalg::Optimize<ptype>::scalar_type;
  struct {
    double k;
    double De;
    scalar operator()(const ptype& psx, ptype& outputs) const {
      size_t n = psx.size();
      assert(n == 1);
      std::vector<scalar> psxk(n);
      std::vector<scalar> output(n);
      psx.get(&(psxk[0]), n, 0);
      auto x = psxk[0];
      auto ex = std::exp(-k * x * 0.5 / std::sqrt(De));
      double fp = De * ((1 - ex) * (1 - ex));
      output[0] = std::sqrt(De) * k * ex * (1 - ex);
      outputs.put(&(output[0]), n, 0);
      return fp;
    }
  } _Morse_residual;

  struct {
    void operator()(ptype& psc, const ptype& psg, std::vector<scalar> shift, bool append = true) const {
      size_t n = psc.size();
      std::vector<scalar> psck(n);
      std::vector<scalar> psgk(n);
      psg.get(&psgk[0], n, 0);
      if (append)
        psc.get(&psck[0], n, 0);
      else
        psck.assign(0, 0);
      for (size_t i = 0; i < n; i++)
        psck[i] -= psgk[i];
      psc.put(&psck[0], n, 0);
    }
  } _Morse_updater;
  bool run(double start, double De) {
    _Morse_residual.k = 1;
    _Morse_residual.De = De;
    constexpr size_t n = 1;
    ptype x(n);
    ptype g(n);
    ptype hg(n);
    const int verbosity = 1;

    molpro::linalg::Optimize<ptype> d("L-BFGS");
    //    IterativeSolver::DIIS<ptype> d;
    d.m_verbosity = verbosity - 1;
    d.m_options["convergence"] = "residual";
    d.m_Wolfe_1 = 1e20;
    d.m_linesearch_tolerance = 0;
    std::vector<scalar> xxx(n, start);
    x.put(&xxx[0], n, 0);
    if (verbosity > 1)
      molpro::cout << "initial guess " << x << std::endl;
    bool converged = false;
    for (int iteration = 1; iteration < 10 && not converged; iteration++) {
      auto value = _Morse_residual(x, g);
      std::vector<scalar> shift;
      shift.push_back(1e-10);
      if (verbosity > 1)
        molpro::cout << "x" << x << std::endl;
      if (verbosity > 1)
        molpro::cout << "g" << g << std::endl;
      if (d.addValue(x, value, g))
        _Morse_updater(x, g, shift);
      converged = d.endIteration(x, g);
      x.get(&xxx[0], n, 0);
      if (verbosity >= 0)
        molpro::cout << "iteration " << iteration << ", error = " << d.errors().front() << ", converged? " << converged
                     << ", value= " << value << std::endl;
    }

    //    x.get(&xxx[0], n, 0);
    //    molpro::cout << "final solution:" ;
    //    for (const auto& xxxx : xxx) molpro::cout << " "<<xxxx; molpro::cout <<std::endl;

    return converged;
  }
};

TEST(Morse_1inf, Optimize) { ASSERT_TRUE(MorseTest().run(1, 1e20)); }
TEST(Morse_110, Optimize) { ASSERT_TRUE(MorseTest().run(1, 10)); }
TEST(Morse_11, Optimize) { ASSERT_TRUE(MorseTest().run(1, 1)); }

class trigTest {

public:
  using scalar = double;
  //using pv = molpro::linalg::PagedArray<scalar>;
  using pv = molpro::linalg::SimpleArray<scalar>;

  size_t n;
  std::string method;
  double hessian;
  scalar initial;

  trigTest(std::string method = "L-BFGS", size_t n = 1, double hessian = 1, double initial = 1)
      : n(n), method(method), hessian(hessian), initial(initial) {}

  scalar residual(const pv& psx, pv& outputs) {
    std::vector<scalar> psxk(n);
    std::vector<scalar> output(n);
    psx.get(psxk.data(), n, 0);
    //    molpro::cout << "residual x: "<<psxk<<std::endl;
    scalar value = 0;
    for (size_t i = 0; i < n; i++) {
      value += 1 - std::cos(scalar(i + 1) * psxk[i]);
      output[i] = scalar(i + 1) * std::sin(scalar(i + 1) * psxk[i]);
    }
    outputs.put(output.data(), n, 0);
    //    molpro::cout << "residual g: "<<output<<std::endl;
    //    molpro::cout << "residual value: "<<value<<std::endl;
    return value;
  }

  void update(pv& psc, const pv& psg) {
    std::vector<scalar> psck(n);
    std::vector<scalar> psgk(n);
    psg.get(psgk.data(), n, 0);
    psc.get(psck.data(), n, 0);
    //    molpro::cout << "update old x: "<<psck<<std::endl;
    //    molpro::cout << "update old g: "<<psgk<<std::endl;
    for (size_t i = 0; i < n; i++)
      psck[i] -= psgk[i] / hessian;
    //    molpro::cout << "update new x: "<<psck<<std::endl;
    psc.put(psck.data(), n, 0);
  }

  bool run() {
    //    std::cout << "optimize with " << method << std::endl;
    molpro::linalg::array::ArrayHandlerIterable<molpro::linalg::SimpleArray<scalar>> handler{};
    molpro::linalg::Optimize<pv> solver(std::regex_replace(method, std::regex("-iterate"), ""));
    solver.m_verbosity = 1;
    solver.m_maxIterations = 50;
    solver.m_thresh = 1e-12;
    //    solver.m_Wolfe_1=.8;
    //    solver.m_linesearch_tolerance = .0001;
    //    std::cout << "Wolfe condition parameters: " << solver.m_Wolfe_1 << ", " << solver.m_Wolfe_2 << std::endl;
    pv g(n);
    pv x(n);
    for (auto i = 0; i < n; i++)
      x.put(&initial, 1, i);
    for (size_t iter = 1; iter <= solver.m_maxIterations; ++iter) {
      auto value = residual(x, g);
      if (solver.m_verbosity > 1)
        molpro::cout << "start iteration " << iter << " value=" << value << "\n x: " << x << "\n g: " << g << std::endl;
      // if (solver.addValue(x, value, g))
      //  update(x, g);
      bool upd = solver.addValue(x, value, g);
      //if (!x.synchronised())
      //  x.sync();
      //if (!g.synchronised())
      //  g.sync();
      if (upd)
        update(x, g);
      // if (solver.endIteration(x, g)) break;
      //      molpro::cout << "before endIteration x "<<x<<std::endl;
      //      molpro::cout << "before endIteration g "<<g<<std::endl;
      upd = solver.endIteration(x, g);
      //      molpro::cout << "after endIteration x "<<x<<std::endl;
      //      molpro::cout << "after endIteration g "<<g<<std::endl;
      //if (!x.synchronised())
      //  x.sync();
      //if (!g.synchronised())
      //  g.sync();
      if (upd)
        break;
    }
    std::cout << "Distance of solution from exact solution: " << std::sqrt(handler.dot(x, x)) << std::endl;
    std::cout << "Error=" << solver.errors().front() << " after " << solver.iterations() << " iterations" << std::endl;
    return std::sqrt(handler.dot(x, x)) < 1e-5 && solver.errors().front() < 1e-5;
  }
};

TEST(Trig_BFGS1, Optimize) { ASSERT_TRUE(trigTest("L-BFGS", 1, 1, 1).run()); }
TEST(Trig_BFGS2, Optimize) { ASSERT_TRUE(trigTest("L-BFGS", 1, 2, 1).run()); }
TEST(Trig_BFGS3, Optimize) { ASSERT_TRUE(trigTest("L-BFGS", 1, 3, 1).run()); }

class optTest {

public:
  using scalar = double;
  //using pv = molpro::linalg::PagedArray<scalar>;
  using pv = molpro::linalg::SimpleArray<scalar>;

protected:
  std::string method;
  std::string name;
  double hessian;
  scalar initial;
  std::vector<scalar> exact;

public:
  optTest(std::string method = "L-BFGS", double hessian = 1, double initial = 1)
      : method(method), hessian(hessian), initial(initial) {}

protected:
  virtual scalar residual(const std::vector<scalar>& psxk, std::vector<scalar>& output) {
    std::size_t n = exact.size();
    scalar value = 0;
    for (size_t i = 0; i < n; i++) {
      value += 1 - std::cos(scalar(i + 1) * psxk[i]);
      output[i] = scalar(i + 1) * std::sin(scalar(i + 1) * psxk[i]);
    }
    return value;
  }
  scalar vresidual(const pv& psx, pv& outputs) {
    std::size_t n = exact.size();
    std::vector<scalar> psxk(n);
    std::vector<scalar> output(n);
    psx.get(psxk.data(), n, 0);
    scalar value = residual(psxk, output);
    outputs.put(output.data(), n, 0);
    return value;
  }

  void update(pv& psc, const pv& psg) {
    std::size_t n = exact.size();
    ASSERT_EQ(n, psc.size());
    std::vector<scalar> psck(n);
    std::vector<scalar> psgk(n);
    psg.get(psgk.data(), n, 0);
    psc.get(psck.data(), n, 0);
    for (size_t i = 0; i < n; i++)
      psck[i] -= psgk[i] / hessian;
    psc.put(psck.data(), n, 0);
  }

public:
  int run(int verbosity = 0) {
    std::size_t n = exact.size();
    if (verbosity > 0)
      std::cout << "optimize " << name << "(" << n << ") with " << method << std::endl;
    molpro::linalg::Optimize<pv> solver(std::regex_replace(method, std::regex("-iterate"), ""));
    solver.m_verbosity = verbosity;
    solver.m_maxIterations = 1000;
    solver.m_thresh = 1e-12;
    //    solver.m_Wolfe_1=.8;
    //    solver.m_linesearch_tolerance = .0001;
    if (verbosity > 0) {
      std::cout << "Wolfe condition parameters: " << solver.m_Wolfe_1 << ", " << solver.m_Wolfe_2 << std::endl;
      std::cout << "initial=" << initial << ", hessian=" << hessian << std::endl;
    }
    pv g(n);
    pv x(n);
    for (auto i = 0; i < n; i++)
      x.put(&initial, 1, i);
    for (size_t iter = 1; iter <= solver.m_maxIterations; ++iter) {
      auto value = vresidual(x, g);
      if (solver.m_verbosity > 1)
        molpro::cout << "start iteration " << iter << " value=" << value << "\n x: " << x << "\n g: " << g << std::endl;
      // if (solver.addValue(x, value, g))
      //  update(x, g);
      bool upd = solver.addValue(x, value, g);
      //if (!x.synchronised())
      //  x.sync();
      //if (!g.synchronised())
      //  g.sync();
      if (upd)
        update(x, g);
      // if (solver.endIteration(x, g)) break;
      upd = solver.endIteration(x, g);
      //if (!x.synchronised())
      //  x.sync();
      //if (!g.synchronised())
      //  g.sync();
      if (upd)
        break;
    }
    std::vector<scalar> xx(n);
    x.get(xx.data(), n, 0);
    scalar dist = 0;
    for (int k = 0; k < n; k++)
      dist += std::pow(xx[k] - exact[k], 2);
    dist = std::sqrt(dist);
    if (verbosity > 0) {
      std::cout << "Distance of solution from exact solution: " << dist << std::endl;
      std::cout << "Error=" << solver.errors().front() << " after " << solver.iterations() << " iterations"
                << std::endl;
    }
    return (dist < 1e-5 && solver.errors().front() < 1e-5) ? solver.iterations() : 1000000;
  }
};

TEST(trigonometric, Optimize) {
  class Test : public optTest {
  public:
    Test(std::string method = "L-BFGS", double hessian = 1, double initial = 1) : optTest(method, hessian, initial) {
      this->name = "trigonometric";
      this->exact.push_back(0);
    }

  private:
    scalar residual(const std::vector<scalar>& psxk, std::vector<scalar>& output) override {
      scalar value = 0;
      size_t n = exact.size();
      for (size_t i = 0; i < n; i++) {
        value += 1 - std::cos(scalar(i + 1) * psxk[i]);
        output[i] = scalar(i + 1) * std::sin(scalar(i + 1) * psxk[i]);
      }
      return value;
    }
  };
  ASSERT_LE(Test().run(0), 6);
}

TEST(Rosenbrock, DISABLED_Optimize) {
  class Test : public optTest {
  public:
    Test(double initial = 2, double hessian = 100) : optTest("L-BFGS", hessian, initial) {
      this->name = "Rosenbrock";
      this->exact.push_back(1);
      this->exact.push_back(1);
    }

  private:
    scalar residual(const std::vector<scalar>& x, std::vector<scalar>& g) override {
      g[0] = -400 * x[0] * (x[1] - std::pow(x[0], 2)) + 2 * (x[0] - 1);
      g[1] = 200 * (x[1] - std::pow(x[0], 2));
      scalar value = 100 * std::pow(x[1] - std::pow(x[0], 2), 2) + std::pow(x[0] - 1, 2);
      //      std::cout << "x=" << x[0] << "," << x[1] << "; g=" << g[0] << "," << g[1] << std::endl;
      //      std::cout << "Rosenbrock residual() at "<<x[0]<<","<<x[1]<<", function="<<value<<",
      //      gradient="<<g[0]<<","<<g[1] << std::endl;
      return value;
    }
  };
  std::map<double, int> expected_iterations; // to catch performance regressions
  expected_iterations[-20] = 63;
  expected_iterations[-2] = 27;
  expected_iterations[-1] = 32;
  expected_iterations[0.01] = 35;
  expected_iterations[1] = 1;
  expected_iterations[2] = 30;
  expected_iterations[3] = 41;
  expected_iterations[8] = 56;
  expected_iterations[30] = 88;
  for (const auto& x : expected_iterations)
    ASSERT_LE(Test(x.first, 800 * x.first * x.first).run(0), x.second);
}
