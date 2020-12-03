#include "create_solver.h"
#include "test.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <numeric>

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
  for (size_t i = 0; i < nX; ++i)
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

std::vector<double> residual(const MatrixXdr& rhs, const std::vector<Rvector>& actions) {
  const auto nroots = actions.size();
  const auto nX = rhs.rows();
  auto errors = std::vector<double>(nroots);
  for (size_t i = 0; i < nroots; ++i) {
    double norm = 0;
    for (size_t j = 0; j < nX; ++j) {
      errors[i] += rhs(i, j) - actions[i][j];
      norm += rhs(i, j) * rhs(i, j);
    }
    if (norm != 0)
      errors[i] /= norm;
    errors[i] = std::sqrt(std::abs(errors[i]));
  }
  return errors;
}

void run_test(const MatrixXdr& mat, const MatrixXdr& rhs, const Update& update, bool augmented_hessian) {
  const auto n_root_max = rhs.cols();
  const auto nX = mat.rows();
  auto reference_solutions = solve_full_problem(mat, rhs, augmented_hessian);
  for (size_t nroot = 1; nroot <= n_root_max; ++nroot) {
    // solve the problem iteratively
    // compare solutions
    std::vector<std::vector<double>> parameters, actions;
    std::vector<int> roots;
    //    for (int root = 0; root < solver->n_roots(); root++) {
    for (int root = 0; root < nroot; root++) {
      parameters.emplace_back(nX);
      actions.emplace_back(nX);
      for (size_t i = 0; i < nX; ++i) {
        parameters[root][i] = reference_solutions[root][i];
        actions[root][i] = rhs(i, root);
      }
      roots.push_back(root);
    }
    //    solver->solution(roots, parameters, actions);
    auto residual_errors = residual(rhs, actions);
    EXPECT_EQ(residual_errors.size(), nroot);
    EXPECT_THAT(residual_errors, ::testing::Each(::testing::Le(1.))); // options->convergence_threshold);
    for (size_t i = 0; i < nroot; ++i) {
      auto norm = std::inner_product(reference_solutions.at(i).begin(), reference_solutions.at(i).end(),
                                     reference_solutions.at(i).begin(), 0.);
      if (norm != 0.) {
        auto overlap_with_reference =
            std::inner_product(parameters.at(i).begin(), parameters.at(i).end(), reference_solutions.at(i).begin(), 0.);
        overlap_with_reference = std::sqrt(std::abs(overlap_with_reference / norm));
        EXPECT_NEAR(overlap_with_reference, 1., 1) << "root = " << i; // options->convergence_threshold.value());
      }
    }
  }
}
} // namespace

TEST(LinearEquations, simple_symmetric_system) {
  double p = 1;
  size_t n_max = 4;
  bool augmented_hessian = false;
  for (size_t n = 1; n < n_max; ++n) {
    auto [mat, rhs] = construct_simple_symmetric_system(n, p);
    auto update = [](const auto& params, const auto& actions, auto n_work) {};
    run_test(mat, rhs, update, augmented_hessian);
  }
}
