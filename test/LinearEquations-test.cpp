#include "create_solver.h"
#include "test.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <numeric>

#include <molpro/linalg/itsolv/CastOptions.h>
#include <molpro/linalg/itsolv/helper.h>

using MatrixXdr = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using Update = std::function<void(const std::vector<Rvector>& params, std::vector<Rvector>& actions, size_t n_work)>;

namespace {
struct LinearEquationsPreconditioner {
  explicit LinearEquationsPreconditioner(const MatrixXdr& mat) {
    const auto n = mat.rows();
    diagonals.resize(n);
    for (size_t i = 0; i < n; ++i)
      diagonals[i] = mat(i, i);
  }

  void operator()(std::vector<Rvector>& residuals, const int nwork) const {
    for (size_t i = 0; i < nwork; ++i) {
      for (size_t j = 0; j < residuals[i].size(); ++j)
        residuals[i][j] /= diagonals[i];
    }
  }

  Rvector diagonals;
};
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
    auto x = Eigen::Map<const Eigen::VectorXd>(params.at(i).data(), params[i].size());
    auto r = Eigen::Map<Eigen::VectorXd>(actions.at(i).data(), params[i].size());
    r = mat * x;
  }
}

auto solve_full_problem(const MatrixXdr& mat, const MatrixXdr& rhs, double augmented_hessian) {
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

std::vector<double> residual(const MatrixXdr& mat, const MatrixXdr& rhs, const std::vector<Rvector>& parameters,
                             std::vector<Rvector>& actions) {
  const auto nroots = parameters.size();
  apply_matrix(mat, parameters, actions, nroots);
  const auto nX = rhs.rows();
  auto errors = std::vector<double>(nroots);
  for (size_t i = 0; i < nroots; ++i) {
    double norm = 0;
    for (size_t j = 0; j < nX; ++j) {
      errors[i] += std::pow(rhs(j, i) - actions[i][j], 2);
      norm += std::pow(rhs(j, i), 2);
    }
    if (norm != 0)
      errors[i] /= norm;
    errors[i] = std::sqrt(std::abs(errors[i]));
  }
  return errors;
}

auto initial_guess(const MatrixXdr& mat, const size_t n_roots) {
  const auto n = mat.rows();
  std::vector<Rvector> g;
  std::vector<Rvector> x;
  std::vector<size_t> guess;
  std::vector<double> diagonals;
  diagonals.reserve(n);
  for (auto i = 0; i < n; i++) {
    diagonals.push_back(mat(i, i));
  }
  for (size_t root = 0; root < n_roots; root++) {
    x.emplace_back(n, 0);
    g.emplace_back(n);
    auto it_min = std::min_element(diagonals.begin(), diagonals.end());
    guess.push_back(std::distance(diagonals.begin(), it_min)); // initial guess
    *std::min_element(diagonals.begin(), diagonals.end()) = 1e99;
    x.back()[guess.back()] = 1; // initial guess
  }
  return std::make_tuple(x, g, guess);
}

auto set_options(std::shared_ptr<ILinearEquations<Rvector, Qvector, Pvector>>& solver, std::shared_ptr<Logger>& logger,
                 const int n, const int nroot, const int np, bool hermitian, double augmented_hessian) {
  auto options = CastOptions::LinearEquations(solver->get_options());
  options->n_roots = nroot;
  options->convergence_threshold = 1.0e-8;
  //    options->norm_thresh = 1.0e-14;
  //    options->svd_thresh = 1.0e-10;
  options->max_size_qspace = std::max(std::min(n, 6 * nroot), std::min(n, std::min(1000, 6 * nroot)) - np);
  options->reset_D = 8;
  options->hermiticity = hermitian;
  options->augmented_hessian = augmented_hessian;
  solver->set_options(options);
  options = CastOptions::LinearEquations(solver->get_options());
  molpro::cout << "convergence threshold = " << options->convergence_threshold.value()
               << ", svd thresh = " << options->svd_thresh.value() << ", norm thresh = " << options->norm_thresh.value()
               << ", max size of Q = " << options->max_size_qspace.value() << ", reset D = " << options->reset_D.value()
               << ", augmented hessian = " << options->augmented_hessian.value() << std::endl;
  logger->max_trace_level = molpro::linalg::itsolv::Logger::Info;
  logger->max_warn_level = molpro::linalg::itsolv::Logger::Error;
  logger->data_dump = true;
  return options;
}

std::vector<std::vector<double>> matrix_to_vector(const MatrixXdr& rhs, const size_t nRHS) {
  auto rhs_vector = std::vector<std::vector<double>>{};
  const auto nX = rhs.rows();
  for (size_t i = 0; i < nRHS; ++i) {
    rhs_vector.emplace_back(nX);
    for (size_t j = 0; j < nX; ++j) {
      rhs_vector.back()[j] = rhs(j, i);
    }
  }
  return rhs_vector;
}
void print_parameters_actions(const std::vector<Rvector>& x, const std::vector<Rvector>& g, size_t nwork) {
  //
  for (size_t i = 0; i < nwork; ++i) {
    std::cout << "parameter_" << i << " = ";
    std::copy(std::begin(x[i]), std::end(x[i]), std::ostream_iterator<double>(std::cout, ", "));
    std::cout << std::endl;
    std::cout << "action_" << i << " = ";
    std::copy(std::begin(g[i]), std::end(g[i]), std::ostream_iterator<double>(std::cout, ", "));
    std::cout << std::endl;
  }
}

void run_test(const MatrixXdr& mat, const MatrixXdr& rhs, const Update& update, double augmented_hessian) {
  auto d = (mat - mat.transpose()).norm();
  bool hermitian = d < 1e-10;
  const auto n_root_max = rhs.cols();
  const auto nX = mat.rows();
  auto reference_solutions = solve_full_problem(mat, rhs, augmented_hessian);
  auto preconditioner = LinearEquationsPreconditioner(mat);
  for (size_t nroot = 1; nroot <= n_root_max; ++nroot) {
    auto rhs_vector = matrix_to_vector(rhs, nroot);
    auto [solver, logger] = molpro::test::create_LinearEquations();
    auto options = set_options(solver, logger, mat.rows(), nroot, 0, hermitian, augmented_hessian);
    solver->add_equations(molpro::linalg::itsolv::cwrap(rhs_vector));
    auto [x, g, guess] = initial_guess(mat, solver->n_roots());
    int nwork = solver->n_roots();
    int max_iter = 100;
    size_t n_iter = 1;
    for (size_t iter = 0; iter < max_iter; ++iter, ++n_iter) {
      apply_matrix(mat, x, g, nwork);
      if (true)
        print_parameters_actions(x, g, nwork);
      nwork = solver->add_vector(x, g);
      if (nwork == 0)
        break;
      solver->report();
      preconditioner(g, nwork);
      nwork = solver->end_iteration(x, g);
      if (nwork == 0)
        break;
    }
    std::vector<std::vector<double>> parameters, actions;
    std::vector<int> roots;
    for (int root = 0; root < solver->n_roots(); root++) {
      //    for (int root = 0; root < nroot; root++) {
      parameters.emplace_back(nX);
      actions.emplace_back(nX);
      if (false)
        for (size_t i = 0; i < nX; ++i) {
          parameters[root][i] = reference_solutions[root][i];
          actions[root][i] = rhs(i, root);
        }
      roots.push_back(root);
    }
    solver->solution(roots, parameters, actions);
    auto residual_errors = residual(mat, rhs, parameters, actions);
    EXPECT_EQ(residual_errors.size(), nroot);
    EXPECT_THAT(residual_errors, ::testing::Each(::testing::Le(options->convergence_threshold.value())));
    for (size_t i = 0; i < nroot; ++i) {
      auto norm = std::inner_product(reference_solutions.at(i).begin(), reference_solutions.at(i).end(),
                                     reference_solutions.at(i).begin(), 0.);
      if (norm != 0.) {
        auto overlap_with_reference =
            std::inner_product(parameters.at(i).begin(), parameters.at(i).end(), reference_solutions.at(i).begin(), 0.);
        overlap_with_reference = std::sqrt(std::abs(overlap_with_reference / norm));
        EXPECT_NEAR(overlap_with_reference, 1., options->convergence_threshold.value()) << "root = " << i;
      }
    }
  }
}
} // namespace

TEST(LinearEquations, simple_symmetric_system) {
  double p = 1;
  size_t n_max = 3;
  double augmented_hessian = 0;
  for (size_t n = 3; n <= n_max; ++n) {
    auto [mat, rhs] = construct_simple_symmetric_system(n, p);
    auto update = [](const auto& params, const auto& actions, auto n_work) {};
    std::cout << "matrix =\n" << mat << std::endl;
    std::cout << "rhs =\n" << rhs << std::endl;
    run_test(mat, rhs, update, augmented_hessian);
  }
}
