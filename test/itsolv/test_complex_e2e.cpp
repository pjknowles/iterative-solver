/*!
 * @file
 * @brief End-to-end check that the solvers run with complex containers.
 *
 * A complex hermitian eigenproblem must converge to real eigenvalues, and must agree with the
 * equivalent real problem when the imaginary parts are switched off -- which is what tells apart a
 * genuinely working complex path from one that merely compiles.
 */
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/itsolv/SolverFactory-implementation.h>

#include <Eigen/Dense>

#include <cmath>
#include <complex>
#include <deque>
#include <map>
#include <string>
#include <vector>

namespace {

template <typename scalar>
struct containers {
  using R = std::vector<scalar>;
  using Q = std::deque<scalar>;
  using P = std::map<size_t, scalar>;
};

/*!
 * @brief A hermitian matrix with a dominant real diagonal.
 *
 * @p mixing scales the off-diagonal block; setting its imaginary part to zero recovers a real
 * symmetric problem, which lets the same problem be solved in both scalar types.
 */
template <typename scalar>
struct HermitianProblem : molpro::linalg::itsolv::Problem<typename containers<scalar>::R> {
  using R = typename containers<scalar>::R;
  size_t n;
  bool complex_coupling;
  HermitianProblem(size_t n, bool complex_coupling) : n(n), complex_coupling(complex_coupling) {}

  scalar matrix(size_t i, size_t j) const {
    if (i == j)
      return scalar(double(i + 1));
    const double re = 0.02 * (double(i) + double(j));
    if constexpr (molpro::linalg::itsolv::is_complex<scalar>{}) {
      const double im = complex_coupling ? 0.02 * (double(j) - double(i)) : 0.0;
      return scalar(re, im);
    } else {
      return scalar(re);
    }
  }

  bool diagonals(R& d) const override {
    for (size_t i = 0; i < d.size(); i++)
      d[i] = matrix(i, i);
    return true;
  }

  void action(const molpro::linalg::itsolv::CVecRef<R>& parameters,
              const molpro::linalg::itsolv::VecRef<R>& actions) const override {
    for (size_t k = 0; k < parameters.size(); k++) {
      const auto& v = parameters[k].get();
      auto& a = actions[k].get();
      for (size_t i = 0; i < a.size(); i++) {
        a[i] = 0;
        for (size_t j = 0; j < a.size(); j++)
          a[i] += matrix(i, j) * v[j];
      }
    }
  }
};

struct Outcome {
  bool converged;
  std::vector<double> eigenvalues;      //!< real parts
  std::vector<double> imaginary_parts;  //!< must vanish for a hermitian problem
  double largest_error;
};

template <typename scalar>
Outcome solve(size_t n, size_t nroot, bool complex_coupling, double threshold, const std::string& options = "") {
  using C = containers<scalar>;
  auto solver = molpro::linalg::itsolv::create_LinearEigensystem<typename C::R, typename C::Q, typename C::P>(
      "Davidson", options);
  solver->set_n_roots(nroot);
  solver->set_convergence_threshold(threshold);
  solver->set_verbosity(0);
  HermitianProblem<scalar> problem(n, complex_coupling);
  std::vector<typename C::R> params(nroot, typename C::R(n, scalar(0))), actions(nroot, typename C::R(n, scalar(0)));
  for (size_t i = 0; i < nroot; ++i)
    params[i][i] = scalar(1);
  Outcome outcome{};
  outcome.converged = solver->solve(params, actions, problem);
  for (size_t i = 0; i < nroot; ++i) {
    const auto e = solver->eigenvalues()[i];
    if constexpr (molpro::linalg::itsolv::is_complex<scalar>{}) {
      outcome.eigenvalues.push_back(e.real());
      outcome.imaginary_parts.push_back(e.imag());
    } else {
      outcome.eigenvalues.push_back(e);
      outcome.imaginary_parts.push_back(0);
    }
  }
  outcome.largest_error = 0;
  for (auto e : solver->errors())
    outcome.largest_error = std::max(outcome.largest_error, double(e));
  return outcome;
}

//! The lowest eigenvalues of the same matrix, computed densely
std::vector<double> reference_eigenvalues(size_t n, size_t nroot, bool complex_coupling) {
  HermitianProblem<std::complex<double>> problem(n, complex_coupling);
  Eigen::MatrixXcd m(n, n);
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j)
      m(i, j) = problem.matrix(i, j);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver(m);
  std::vector<double> values;
  for (size_t i = 0; i < nroot; ++i)
    values.push_back(solver.eigenvalues()(i));
  return values;
}

constexpr size_t n = 14;
constexpr size_t nroot = 3;

} // namespace

template class molpro::linalg::itsolv::SolverFactory<std::vector<std::complex<double>>,
                                                     std::deque<std::complex<double>>,
                                                     std::map<size_t, std::complex<double>>>;

//! A complex hermitian eigenproblem converges to the eigenvalues of the dense matrix, and they are real
TEST(complex_e2e, davidson_hermitian) {
  const auto outcome = solve<std::complex<double>>(n, nroot, true, 1e-11);
  EXPECT_TRUE(outcome.converged);
  EXPECT_LT(outcome.largest_error, 1e-11);
  const auto reference = reference_eigenvalues(n, nroot, true);
  for (size_t i = 0; i < nroot; ++i) {
    EXPECT_NEAR(outcome.eigenvalues[i], reference[i], 1e-10) << "root " << i;
    EXPECT_NEAR(outcome.imaginary_parts[i], 0.0, 1e-12) << "root " << i << " of a hermitian problem must be real";
  }
}

//! With the imaginary coupling switched off, the complex run must reproduce the real one exactly
TEST(complex_e2e, davidson_matches_real_when_coupling_is_real) {
  const auto in_complex = solve<std::complex<double>>(n, nroot, false, 1e-11);
  const auto in_double = solve<double>(n, nroot, false, 1e-11);
  ASSERT_TRUE(in_complex.converged);
  ASSERT_TRUE(in_double.converged);
  for (size_t i = 0; i < nroot; ++i) {
    EXPECT_NEAR(in_complex.eigenvalues[i], in_double.eigenvalues[i], 1e-10) << "root " << i;
    EXPECT_NEAR(in_complex.imaginary_parts[i], 0.0, 1e-12) << "root " << i;
  }
}

//! The imaginary coupling changes the spectrum, so the test above is not vacuous
TEST(complex_e2e, imaginary_coupling_changes_the_spectrum) {
  const auto with_imaginary = reference_eigenvalues(n, nroot, true);
  const auto without = reference_eigenvalues(n, nroot, false);
  bool differ = false;
  for (size_t i = 0; i < nroot; ++i)
    differ = differ || std::abs(with_imaginary[i] - without[i]) > 1e-6;
  EXPECT_TRUE(differ) << "the two problems must not be the same, or the comparison above proves nothing";
}

namespace {
//! Checks a complex hermitian Davidson run with the given options against the dense eigenvalues
void check_davidson_hermitian(const std::string& options) {
  const auto outcome = solve<std::complex<double>>(n, nroot, true, 1e-10, options);
  EXPECT_TRUE(outcome.converged);
  EXPECT_LT(outcome.largest_error, 1e-10);
  const auto reference = reference_eigenvalues(n, nroot, true);
  for (size_t i = 0; i < nroot; ++i) {
    EXPECT_NEAR(outcome.eigenvalues[i], reference[i], 1e-9) << "root " << i;
    EXPECT_NEAR(outcome.imaginary_parts[i], 0.0, 1e-12) << "root " << i << " of a hermitian problem must be real";
  }
}
} // namespace

//! A Q space too small to hold the expansion forces construction of a D space from deleted Q vectors
TEST(complex_e2e, davidson_hermitian_dspace) { check_davidson_hermitian("max_size_qspace=4"); }

//! As above, with the D space also being reset periodically
TEST(complex_e2e, davidson_hermitian_dspace_reset) { check_davidson_hermitian("max_size_qspace=4,reset_D=3"); }

namespace {

using complex_scalar = std::complex<double>;
using CC = containers<complex_scalar>;

/*!
 * @brief A hermitian positive definite operator with a known right-hand side and a known minimum.
 *
 * Serves the three solver families that the eigensolver tests above do not reach: linear equations
 * A x = b, the non-linear equations A(x - 1) = 0, and the minimisation of the real quadratic form
 * 1/2 <x-1|A|x-1>. Each has a solution that can be checked without trusting the solver.
 */
struct HermitianPositiveDefinite : molpro::linalg::itsolv::Problem<CC::R> {
  size_t n;
  explicit HermitianPositiveDefinite(size_t n) : n(n) {}

  complex_scalar matrix(size_t i, size_t j) const {
    if (i == j)
      return complex_scalar(double(i + 1), 0);
    return complex_scalar(0.02 * (double(i) + double(j)), 0.02 * (double(j) - double(i)));
  }

  //! The right-hand side, chosen so that its phase varies along the vector
  complex_scalar rhs(size_t i) const { return complex_scalar(0.1 * double(i + 1), -0.05 * double(i)); }

  Eigen::MatrixXcd dense() const {
    Eigen::MatrixXcd m(n, n);
    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < n; ++j)
        m(i, j) = matrix(i, j);
    return m;
  }

  bool diagonals(CC::R& d) const override {
    for (size_t i = 0; i < d.size(); i++)
      d[i] = matrix(i, i);
    return true;
  }

  void action(const molpro::linalg::itsolv::CVecRef<CC::R>& parameters,
              const molpro::linalg::itsolv::VecRef<CC::R>& actions) const override {
    for (size_t k = 0; k < parameters.size(); k++) {
      const auto& v = parameters[k].get();
      auto& a = actions[k].get();
      for (size_t i = 0; i < a.size(); i++) {
        a[i] = 0;
        for (size_t j = 0; j < a.size(); j++)
          a[i] += matrix(i, j) * v[j];
      }
    }
  }

  bool RHS(CC::R& b, unsigned int instance) const override {
    if (instance > 0)
      return false;
    for (size_t i = 0; i < b.size(); ++i)
      b[i] = rhs(i);
    return true;
  }

  //! Residual A(x - 1), whose zero is x = 1; the value is the real quadratic form 1/2 <x-1|A|x-1>
  complex_scalar residual(const CC::R& x, CC::R& a) const override {
    double value = 0;
    for (size_t i = 0; i < a.size(); i++) {
      a[i] = 0;
      for (size_t j = 0; j < a.size(); j++)
        a[i] += matrix(i, j) * (x[j] - complex_scalar(1));
      value += 0.5 * (std::conj(x[i] - complex_scalar(1)) * a[i]).real();
    }
    // the objective of a minimisation is real; the complex scalar type carries it with no imaginary part
    return complex_scalar(value, 0);
  }
};

//! Largest deviation of x from the vector of ones
double distance_from_ones(const CC::R& x) {
  double worst = 0;
  for (const auto& e : x)
    worst = std::max(worst, std::abs(e - complex_scalar(1)));
  return worst;
}

} // namespace

//! Complex linear equations: the solution must match a dense solve of the same hermitian system
TEST(complex_e2e, linear_equations) {
  constexpr size_t dim = 12;
  auto solver = molpro::linalg::itsolv::create_LinearEquations<CC::R, CC::Q, CC::P>("Davidson");
  solver->set_n_roots(1);
  solver->set_convergence_threshold(1e-11);
  solver->set_verbosity(0);
  HermitianPositiveDefinite problem(dim);
  std::vector<CC::R> x(1, CC::R(dim, complex_scalar(0))), g(1, CC::R(dim, complex_scalar(0)));
  // an explicit starting vector: generate_initial_guess orders the diagonal, which select() refuses to
  // do for a complex element type
  x[0][0] = complex_scalar(1);
  EXPECT_TRUE(solver->solve(x, g, problem));
  solver->solution(std::vector<int>{0}, x, g);

  Eigen::VectorXcd b(dim);
  for (size_t i = 0; i < dim; ++i)
    b(i) = problem.rhs(i);
  const Eigen::VectorXcd reference = problem.dense().ldlt().solve(b);
  double worst = 0;
  for (size_t i = 0; i < dim; ++i)
    worst = std::max(worst, std::abs(x[0][i] - reference(i)));
  EXPECT_LT(worst, 1e-9) << "solution differs from the dense solve";
}

//! Complex non-linear equations by DIIS: the residual A(x - 1) has its zero at x = 1
TEST(complex_e2e, diis_nonlinear_equations) {
  constexpr size_t dim = 10;
  auto solver = molpro::linalg::itsolv::create_NonLinearEquations<CC::R, CC::Q, CC::P>("DIIS");
  solver->set_convergence_threshold(1e-11);
  solver->set_verbosity(0);
  HermitianPositiveDefinite problem(dim);
  CC::R x(dim, complex_scalar(0)), g(dim, complex_scalar(0));
  EXPECT_TRUE(solver->solve(x, g, problem));
  EXPECT_LT(distance_from_ones(x), 1e-9);
}

/*!
 * @brief Complex minimisation by BFGS, of a real objective over complex parameters.
 *
 * The line search orders function values, so only the real part of the value is used; that is
 * meaningful here because the quadratic form of a hermitian operator is real.
 */
TEST(complex_e2e, bfgs_real_objective_of_complex_parameters) {
  constexpr size_t dim = 10;
  auto solver = molpro::linalg::itsolv::create_Optimize<CC::R, CC::Q, CC::P>("BFGS");
  solver->set_convergence_threshold(1e-8);
  solver->set_max_iter(200);
  solver->set_verbosity(0);
  HermitianPositiveDefinite problem(dim);
  CC::R x(dim, complex_scalar(0)), g(dim, complex_scalar(0));
  solver->solve(x, g, problem);
  EXPECT_LT(distance_from_ones(x), 1e-6);
}
