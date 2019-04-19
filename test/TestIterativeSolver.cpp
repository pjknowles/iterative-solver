#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "test.h"

#include <Eigen/Dense>
#include <limits>
#include <cmath>
#include <numeric>
#include <vector>
#include "IterativeSolver.h"
#include "SimpleVector.h"
#include "PagedVector.h"

TEST(TestIterativeSolver, small_eigenproblem) {
  for (size_t n = 1; n < 20; n++) {
    for (size_t nroot = 1; nroot <= n && nroot < 10; nroot++) {
      Eigen::MatrixXd m(n, n);
      for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
          m(i, j) = 1 + (i + j) * std::sqrt(double(i + j));
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> denseSolver(m);
      auto val = denseSolver.eigenvalues();

      LinearAlgebra::SimpleVector<double> mm(n);
      std::vector<LinearAlgebra::SimpleVector<double> > x, g;
      IterativeSolver::LinearEigensystem<LinearAlgebra::SimpleVector<double> > solver;
      solver.m_verbosity = -1;
      solver.setThresholds(1e-13);
      if (solver.m_verbosity > 0) std::cout << "Test n=" << n << ", nroot=" << nroot << std::endl;
      for (size_t root = 0; root < nroot; root++) {
        x.emplace_back(n);
        x.back().scal(0);
        x.back()[root] = 1;
        g.emplace_back(n);
      }
      for (size_t iter = 0; iter < n+1; iter++) {
        for (size_t root = 0; root < x.size(); root++) {
          g[root].scal(0);
          if (solver.active()[root])
            for (size_t i = 0; i < n; i++)
              for (size_t j = 0; j < n; j++)
                g[root][j] += m(j, i) * x[root][i];
        }

//        std::cout << "eigenvector "<<0<<active[0]<<" before addVector"; for (size_t i = 0; i < n; i++) std::cout << " " << x[0][i]; std::cout << std::endl;
        solver.addVector(x, g);
        for (size_t root = 0; root < x.size(); root++) {
          if (solver.m_verbosity > 1) {
            std::cout << "eigenvector "<<root<<" before update";
            for (size_t i = 0; i < n; i++) std::cout << " " << x[root][i];
            std::cout << std::endl;
          }
          if (solver.active()[root]) {
            for (size_t i = 0; i < n; i++)
              x[root][i] -= g[root][i] / (m(i, i) - solver.eigenvalues()[root] + 1e-13);
            if (solver.m_verbosity > 2) {
              std::cout << "residual "<<root<<" ";
              for (size_t i = 0; i < n; i++) std::cout << " " << g[root][i];
              std::cout << std::endl;
            }
            if (solver.m_verbosity > 1) {
              std::cout << "eigenvector "<<root<<" ";
              for (size_t i = 0; i < n; i++) std::cout << " " << x[root][i];
              std::cout << std::endl;
            }
          }
        }
//        std::cout << "eigenvector "<<0<<active[0]<<" before endIteration"; for (size_t i = 0; i < n; i++) std::cout << " " << x[0][i]; std::cout << std::endl;
//        auto conv = (solver.endIteration(x, g, active));
//        std::cout << "eigenvector "<<0<<active[0]<<" after endIteration"; for (size_t i = 0; i < n; i++) std::cout << " " << x[0][i]; std::cout << std::endl;
//        if (conv) break;
        if (solver.endIteration(x, g)) break;
      }
//  std::cout << "Error={ "; for (const auto& e : solver.errors()) std::cout << e << " "; std::cout << "} after " << solver.iterations() << " iterations" << std::endl;
//  std::cout << "Actual eigenvalues\n"<<val<<std::endl;
//      EXPECT_THAT(active,
//                  ::testing::Pointwise(::testing::Eq(), std::vector<bool>(nroot, false)));
      EXPECT_THAT(solver.errors(),
                  ::testing::Pointwise(::testing::DoubleNear(1e-10), std::vector<double>(nroot, double(0))));
      EXPECT_THAT(solver.eigenvalues(),
                  ::testing::Pointwise(::testing::DoubleNear(1e-10),
                                       std::vector<double>(val.data(), val.data() + nroot)));
      for (size_t root = 0; root < solver.m_roots; root++) {
        if (solver.m_verbosity > 1) {
          std::cout << "eigenvector "<<root<<" active="<<solver.active()[root]<<" converged="<<solver.errors()[root]<<":";
          for (size_t i = 0; i < n; i++) std::cout << " " << x[root][i];
          std::cout << std::endl;
        }
        std::vector<double> r(n);
        g[root].get(r.data(), n, 0);
        EXPECT_THAT(r, ::testing::Pointwise(::testing::DoubleNear(1e-5), std::vector<double>(n, double(0))));
        if (solver.m_verbosity>1)
          for (size_t soot = 0; soot <= root; soot++)
            std::cout << "Eigenvector overlap "<<root<<" "<<soot<<" "<<x[root].dot(x[soot])<<std::endl;
        for (size_t soot = 0; soot < root; soot++)
          EXPECT_LE (std::abs(x[root].dot(x[soot])), 1e-8); // can't expect exact orthogonality when last thing might have been an update
        EXPECT_THAT (std::abs(x[root].dot(x[root])), ::testing::DoubleNear(1,1e-10));
      }
    }

  }
}

TEST(TestIterativeSolver, linear_equations) {
  for (size_t n = 1; n < 20; n++) {
    for (size_t nroot = 1; nroot <= n && nroot < 10; nroot++) {
      Eigen::MatrixXd m(n, n);
      for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
          m(i, j) = 1 + (i + j) * std::sqrt(double(i + j));
      for (size_t i = 0; i < n; i++)
        m(i,i) += i;

      LinearAlgebra::SimpleVector<double> mm(n);
      std::vector<LinearAlgebra::SimpleVector<double> > x, g, rhs;
      for (size_t root = 0; root < nroot; root++) {
        x.emplace_back(n);
        x.back().scal(0);
        x.back()[root] = 1;
        g.emplace_back(n);
        rhs.emplace_back(n);
        rhs.back().scal(0);
        rhs.back()[root] = 1;
        Eigen::VectorXd erhs(n);
        rhs[root].get(&erhs(0),n,0);
        auto trueSolution = m.colPivHouseholderQr().solve(erhs).eval();
        rhs.back()[root]=1/trueSolution(root);
      }
      IterativeSolver::LinearEquations<LinearAlgebra::SimpleVector<double> > solver(rhs);
      solver.m_verbosity = 0;
      solver.setThresholds(1e-13);
      if (solver.m_verbosity > 0) std::cout << "Test n=" << n << ", nroot=" << nroot << std::endl;
      if (solver.m_verbosity > 1) std::cout << "Matrix:\n"<<m<<std::endl;
      for (size_t iter = 0; iter < n+1; iter++) {
        for (size_t root = 0; root < x.size(); root++) {
          g[root].scal(0);
          if (solver.active()[root])
            for (size_t i = 0; i < n; i++)
              for (size_t j = 0; j < n; j++)
                g[root][j] += m(j, i) * x[root][i];
        }

//        std::cout << "solution "<<0<<solver.active()[0]<<" before addVector"; for (size_t i = 0; i < n; i++) std::cout << " " << x[0][i]; std::cout << std::endl;
        solver.addVector(x, g);
        for (size_t root = 0; root < x.size(); root++) {
          if (solver.m_verbosity > 1) {
            std::cout << "solution "<<root<<" before update";
            for (size_t i = 0; i < n; i++) std::cout << " " << x[root][i];
            std::cout << std::endl;
          }
          if (solver.active()[root]) {
            for (size_t i = 0; i < n; i++)
              x[root][i] -= g[root][i] /
                  (m(i, i)
//                  - solver.eigenvalues()[root]
                  + 1e-13);
            if (solver.m_verbosity > 2) {
              std::cout << "residual "<<root<<" ";
              for (size_t i = 0; i < n; i++) std::cout << " " << g[root][i];
              std::cout << std::endl;
            }
            if (solver.m_verbosity > 1) {
              std::cout << "solution "<<root<<" ";
              for (size_t i = 0; i < n; i++) std::cout << " " << x[root][i];
              std::cout << std::endl;
            }
          }
        }
//        std::cout << "eigenvector "<<0<<active[0]<<" before endIteration"; for (size_t i = 0; i < n; i++) std::cout << " " << x[0][i]; std::cout << std::endl;
//        auto conv = (solver.endIteration(x, g, active));
//        std::cout << "eigenvector "<<0<<active[0]<<" after endIteration"; for (size_t i = 0; i < n; i++) std::cout << " " << x[0][i]; std::cout << std::endl;
//        if (conv) break;
        if (solver.endIteration(x, g)) break;
      }
//      EXPECT_THAT(solver.active(),
//                  ::testing::Pointwise(::testing::Eq(), std::vector<bool>(nroot, false)));
      EXPECT_THAT(solver.errors(),
                  ::testing::Pointwise(::testing::DoubleNear(1e-10), std::vector<double>(nroot, double(0))));
//                                       std::vector<double>(val.data(), val.data() + nroot)));
      for (size_t root = 0; root < solver.m_roots; root++) {
        if (solver.m_verbosity > 1) {
          std::cout << "solution "<<root<<" active="<<solver.active()[root]<<" converged="<<solver.errors()[root]<<":";
          for (size_t i = 0; i < n; i++) std::cout << " " << x[root][i];
          std::cout << std::endl;
        }
        Eigen::VectorXd erhs(n);
        rhs[root].get(&erhs(0),n,0);
        auto trueSolution = m.colPivHouseholderQr().solve(erhs).eval();
//        std::cout << "RHS\n"<<erhs<<std::endl;
//        std::cout << "trueSolution\n"<<trueSolution<<std::endl;
        std::vector<double> r(n);
        x[root].get(r.data(), n, 0);
        EXPECT_THAT(r, ::testing::Pointwise(::testing::DoubleNear(1e-5), std::vector<double>(&trueSolution(0),&trueSolution(0)+n)));
      }
    }

  }
}

#include <regex>
class RosenbrockTest {
 public:
  using ptype = LinearAlgebra::SimpleVector<double>;
  using scalar = typename IterativeSolver::Optimize<ptype>::scalar_type;
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
      return (Rosenbrock_a - psxk[0]) * (Rosenbrock_a - psxk[0]) + Rosenbrock_b * (psxk[1] - psxk[0] * psxk[0]) * (psxk[1] - psxk[0] * psxk[0]);
    }
  } _Rosenbrock_residual;

  struct {
    bool approximateHessian = true;
    void operator()(ptype& psc,
                    const ptype& psg,
                    std::vector<scalar> shift,
                    bool append = true) const {
      size_t n = 2;
      std::vector<scalar> psck(n);
      std::vector<scalar> psgk(n);
      psg.get(&psgk[0], n, 0);
      if (append)
        psc.get(&psck[0], n, 0);
      else
        psck.assign(0, 0);
//      xout << "Rosenbrock updater, initial psc=" << psc << std::endl;
//      xout << "Rosenbrock updater, initial psg=" << psg << std::endl;
//      xout << "Rosenbrock updater, approximateHessian=" << approximateHessian << std::endl;
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
        dx =
            -psgk[0] / (4 + 8 * Rosenbrock_b * (x * x - y))
                - psgk[1] / (4 + 4 * Rosenbrock_b * (x * x - y));
        dy =
            -psgk[1] / (4 * Rosenbrock_b)
                - psgk[1] * x * x / (1 + 2 * Rosenbrock_b * (x * x - y))
                - psgk[0] / (4 + 4 * Rosenbrock_b * (x * x - y));
      }
      psck[0] += dx;
      psck[1] += dy;
      psc.put(&psck[0], n, 0);
//      xout << "Rosenbrock updater, new psc=" << psc << std::endl;
    }
  } _Rosenbrock_updater;
  bool run(std::string method) {

    ptype x(2);
    ptype g(2);
    ptype hg(2);
    const double difficulty = .9;
    const int verbosity = 1;

    if (verbosity >= 0) xout << "Test Optimize, method=" << method << ", difficulty=" << difficulty << std::endl;
    IterativeSolver::Optimize<ptype> d(regex_replace(method, std::regex("-.*"), ""));
    d.m_verbosity = verbosity - 1;
    d.m_options["convergence"] = "residual";
    std::vector<scalar> xxx(2);
    xxx[0] = xxx[1] = 1 - difficulty; // initial guess
    x.put(&xxx[0], 2, 0);
//    xout << "initial guess " << x << std::endl;
    _Rosenbrock_updater.approximateHessian = true;
    bool converged = false;
    for (int iteration = 1; iteration < 50 && not converged; iteration++) {
//      xout << "start of iteration " << iteration << std::endl;
      auto value = _Rosenbrock_residual(x, g);
//      xout << "x: " << x;
//      xout << "g: " << g;
      std::vector<scalar> shift;
      shift.push_back(1e-10);
      if (method == "null-iterate") {
        hg.scal(0);
        _Rosenbrock_updater(hg, g, shift);
//        xout << "hg: " << hg;
        converged = d.iterate(x, g, hg, value);
      } else {
        d.addVector(x, g);
//        xout << "x before update " << x << std::endl;
//        xout << "g before update " << g << std::endl;
        _Rosenbrock_updater(x, g, shift);
//        xout << "x after update " << x << std::endl;
        converged = d.endIteration(x, g);
      }
      x.get(&xxx[0], 2, 0);
//      if (verbosity > 2)
//        xout << "new x after iterate " << x << std::endl;
      if (verbosity >= 0)
        xout << "iteration " << iteration
             << ", Distance from solution = " << std::sqrt(
            (xxx[0] - Rosenbrock_a) * (xxx[0] - Rosenbrock_a) + (xxx[1] - Rosenbrock_a) * (xxx[1] - Rosenbrock_a))
             << ", error = " << d.errors().front()
             << ", converged? " << converged
             << ", value= " << value
             << std::endl;
//   xout <<"end of iteration "<<iteration<<std::endl;
    }

    x.get(&xxx[0], 2, 0);
    xout << "Distance from solution = " << std::sqrt(
        (xxx[0] - Rosenbrock_a) * (xxx[0] - Rosenbrock_a) + (xxx[1] - Rosenbrock_a) * (xxx[1] - Rosenbrock_a))
         << std::endl;

    return converged;
  }
};

TEST(Rosenbrock_BFGS,Optimize) { ASSERT_TRUE(RosenbrockTest().run("BFGS")); }
TEST(Rosenbrock_null,Optimize) { ASSERT_TRUE(RosenbrockTest().run("null")); }
