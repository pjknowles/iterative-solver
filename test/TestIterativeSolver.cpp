#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <limits>
#include <cmath>
#include <numeric>
#include <vector>
#include "IterativeSolver.h"
#include "SimpleVector.h"

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


