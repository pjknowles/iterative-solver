#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <limits>
#include <cmath>
#include <numeric>
#include <vector>
#include "IterativeSolver.h"
#include "PagedVector.h"

TEST(TestIterativeSolver, two_by_two) {
  constexpr size_t n = 4, nroot=1;
  Eigen::MatrixXd m(n,n);
  for (auto i=0; i<n; i++)
    for (auto j=0; j<n; j++)
      m(i,j) = 1+i+j;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> denseSolver(m);
  auto val = denseSolver.eigenvalues();
  auto vec = denseSolver.eigenvectors();
  auto threshold = 20 * std::numeric_limits<double>::epsilon();
  auto overlap = vec.transpose() * vec;
  auto o0 = overlap;
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < i; j++)
      ASSERT_LE (std::fabs(o0(i, j)), threshold);
  }

  LinearAlgebra::PagedVector<double> mm(n);
  std::vector<LinearAlgebra::PagedVector<double> > x,g;
  LinearAlgebra::LinearEigensystem<LinearAlgebra::PagedVector<double> > solver;
  solver.m_verbosity=2;
  std::vector<bool> active;
  for (auto root=0; root<nroot; root++) {
    x.emplace_back(n); x.back().scal(0); x.back()[root]=1; g.emplace_back(n); active.push_back(true);
  }
  for (auto iter = 0; iter < n; iter++) {
    for (auto root=0; root<x.size(); root++) {
      g[root].scal(0);
      for (auto i=0; i<n; i++)
        for (auto j=0; j<n; j++)
          g[root][j] += m(j,i) * x[root][i];
    }

    solver.addVector(x, g, active);
      for (auto root=0; root<x.size(); root++)
        for (auto i=0; i<n; i++)
          for (auto j=0; j<n; j++)
            x[root][j] +=  x[root][i] - g[root][i]/(solver.eigenvalues()[root]-m(i,i)+1e-13);
    if (solver.endIteration(x, g, active)) break;
  }
  std::cout << "Error={ ";
  for (const auto& e : solver.errors()) std::cout << e << " ";
  std::cout << "} after " << solver.iterations() << " iterations" << std::endl;
  for (size_t root = 0; root < solver.m_roots; root++) {
    std::cout << "Eigenvector: (norm=" << std::sqrt(x[root].dot(x[root])) << "): ";
    for (size_t k = 0; k < n; k++) std::cout << " " << (x[root])[k];
    std::cout << std::endl;
  }
}



