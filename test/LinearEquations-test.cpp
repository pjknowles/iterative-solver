#include "create_solver.h"
#include "test.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <molpro/linalg/itsolv/helper.h>

using MatrixXdr = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using Update = std::function<void(const std::vector<Rvector>& params, std::vector<Rvector>& actions, size_t n_work)>;

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

auto solve_full_problem(const MatrixXdr& mat, const MatrixXdr& rhs, bool augmented_hessian) {
  const auto nX = mat.rows(), nroot = rhs.cols();
  std::vector<double> solution, eigenvalues, matrix, metric, rhs_flat;
  matrix.resize(nX * nX);
  metric.resize(nX * nX);
  rhs_flat.resize(nX * nroot);
  for (size_t i = 0, ij = 0; i < nX; ++i)
    metric[i * (nX + 1)] = 1;
  for (size_t i = 0, ij = 0; i < nX; ++i)
    for (size_t j = 0; j < nX; ++j, ++ij)
      matrix[ij] = mat(i, j);
  for (size_t i = 0, ij = 0; i < nX; ++i)
    for (size_t j = 0; j < nroot; ++j, ++ij)
      rhs_flat[ij] = rhs(i, j);
  molpro::linalg::itsolv::solve_LinearEquations(solution, eigenvalues, matrix, metric, rhs_flat, mat.rows(), rhs.cols(),
                                                augmented_hessian, 1.0e-14, 0);
  std::vector<std::vector<double>> solutions_mat;
  for (size_t i = 0, ij = 0; i < nroot; ++i) {
    solutions_mat.emplace_back(nX);
    for (size_t j = 0; j < nX; ++j, ++ij) {
      solutions_mat.back()[j] = solution[ij];
    }
  }
  return solutions_mat;
}

void run_test(const MatrixXdr& mat, const MatrixXdr& rhs, const Update& update, bool augmented_hessian) {
  const auto n_root_max = rhs.cols();
  auto reference_solutions = solve_full_problem(mat, rhs, augmented_hessian);
  for (size_t nroot = 1; nroot <= n_root_max; ++nroot) {
    // solve the problem iteratively
    // compare solutions
  }
}
} // namespace

TEST(LinearEquations, simple_symmetric_system) {
  double p = 1;
  size_t n_max = 4;
  for (size_t n = 0; n < n_max; ++n) {
    auto [mat, rhs] = construct_simple_symmetric_system(n, p);
    auto update = [](const auto& params, const auto& actions, auto n_work) {};
    // choose matrix apply function
    // choose update functon
    // run test
    run_test(mat, rhs, update, false);
  }
}
