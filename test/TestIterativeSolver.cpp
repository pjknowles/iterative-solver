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
  for (size_t n = 5; n < 6; n++) {
    for (size_t nroot = 3; nroot <= n; nroot++) {
      Eigen::MatrixXd m(n, n);
      for (auto i = 0; i < n; i++)
        for (auto j = 0; j < n; j++)
          m(i, j) = 1 + (i + j) * (i + j);
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> denseSolver(m);
      auto val = denseSolver.eigenvalues();

      LinearAlgebra::SimpleVector<double> mm(n);
      std::vector<LinearAlgebra::SimpleVector<double> > x, g;
      LinearAlgebra::LinearEigensystem<LinearAlgebra::SimpleVector<double> > solver;
      solver.m_verbosity = 2;
      if (solver.m_verbosity > 0) std::cout << "Test n=" << n << ", nroot=" << nroot << std::endl;
      std::vector<bool> active;
      for (auto root = 0; root < nroot; root++) {
        x.emplace_back(n);
        x.back().scal(0);
        x.back()[root] = 1;
        g.emplace_back(n);
        active.push_back(true);
      }
      for (auto iter = 0; iter < n; iter++) {
        for (auto root = 0; root < x.size(); root++) {
          g[root].scal(0);
          if (active[root])
            for (auto i = 0; i < n; i++)
              for (auto j = 0; j < n; j++)
                g[root][j] += m(j, i) * x[root][i];
        }

        solver.addVector(x, g, active);
        for (auto root = 0; root < x.size(); root++) {
          if (solver.m_verbosity > 1) {
            std::cout << "eigenvector before update";
            for (auto i = 0; i < n; i++) std::cout << " " << x[root][i];
            std::cout << std::endl;
          }
          if (active[root]) {
            for (auto i = 0; i < n; i++)
              x[root][i] -= g[root][i] / (m(i, i) - solver.eigenvalues()[root] + 1e-13);
            if (solver.m_verbosity > 2) {
              std::cout << "residual";
              for (auto i = 0; i < n; i++) std::cout << " " << g[root][i];
              std::cout << std::endl;
            }
            if (solver.m_verbosity > 1) {
              std::cout << "eigenvector";
              for (auto i = 0; i < n; i++) std::cout << " " << x[root][i];
              std::cout << std::endl;
            }
          }
        }
        if (solver.endIteration(x, g, active)) break;
      }
//  std::cout << "Error={ "; for (const auto& e : solver.errors()) std::cout << e << " "; std::cout << "} after " << solver.iterations() << " iterations" << std::endl;
//  std::cout << "Actual eigenvalues\n"<<val<<std::endl;
      EXPECT_THAT(solver.errors(),
                  ::testing::Pointwise(::testing::DoubleNear(1e-12), std::vector<double>(nroot, double(0))));
      EXPECT_THAT(solver.eigenvalues(),
                  ::testing::Pointwise(::testing::DoubleNear(1e-12),
                                       std::vector<double>(val.data(), val.data() + nroot)));
      for (size_t root = 0; root < solver.m_roots; root++) {
        std::vector<double> r(n);
        g[root].get(r.data(), n, 0);
        EXPECT_THAT(r, ::testing::Pointwise(::testing::DoubleNear(1e-12), std::vector<double>(n, double(0))));
        if (solver.m_verbosity>1)
        for (size_t soot = 0; soot <= root; soot++)
          std::cout << "Eigenvector overlap "<<root<<" "<<soot<<" "<<x[root].dot(x[soot])<<std::endl;
        for (size_t soot = 0; soot < root; soot++)
          EXPECT_LE (std::abs(x[root].dot(x[soot])), 1e-12);
        EXPECT_THAT (std::abs(x[root].dot(x[root])), ::testing::DoubleNear(1,1e-12));
      }
    }

  }
}


