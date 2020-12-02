#include "create_solver.h"
#include "test.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <molpro/linalg/itsolv/helper.h>

using MatrixXdr = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

namespace {
/*
 * Simple symmetric system of linear equations
 * A.x = b
 * A_ij = If[ i==j, i] + p
 * b_ij = i + j
 *
 * Solutions are:
 * x_ij = If[ i == 0, 1 - n + j / p, 1]
 * ^ solutions stored in columns. Consistent with Mathematica and Eigen, but we store them as rows.
 */
auto construct_simple_symmetric_system(size_t n, double p) {
  auto mat = MatrixXdr(n, n);
  auto rhs = MatrixXdr(n, n);
  mat.fill(p);
  for (size_t i = 0; i < n; ++i) {
    mat(i, i) += i;
    for (size_t j = 0; j < n; ++j) {
      rhs(i, j) = i + j;
    }
  }
  return std::make_tuple(mat, rhs);
}

void apply_matrix(const MatrixXdr& mat, const std::vector<Rvector>& params, std::vector<Rvector>& actions,
                  size_t n_work) {
  for (size_t i = 0; i < n_work; ++i) {
    auto x = Eigen::Map<const Eigen::VectorXd>(params.at(i).data(), 1, params[i].size());
    auto r = Eigen::Map<Eigen::VectorXd>(actions.at(i).data(), 1, params[i].size());
    r += mat * x;
  }
}
} // namespace

TEST(LinearEquations, simple_symmetric_system) {
  double p = 1;
  size_t n_max = 4;
  for (size_t n = 0; n < n_max; ++n) {
    auto [mat, rhs] = construct_simple_symmetric_system(n, p);
  }
}
