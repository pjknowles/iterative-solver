/*!
 * @file
 * @brief End-to-end check that the solvers themselves run in a precision other than double.
 *
 * The same eigenproblem is solved with std::vector<double> and std::vector<long double> containers.
 * Only the extended-precision run can drive the residuals below the double-precision floor, which is
 * what makes this a test of the arithmetic rather than of the interface.
 */
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/itsolv/SolverFactory-implementation.h>

#include <cmath>
#include <limits>
#include <deque>
#include <map>
#include <vector>

namespace {

/*!
 * @brief Whether long double is genuinely wider than double on this platform.
 *
 * It is not everywhere: on arm64 macOS, for instance, both are IEEE binary64, and a run in long
 * double can then no more reach below the double-precision floor than a run in double can.
 */
constexpr bool long_double_is_extended =
    std::numeric_limits<long double>::epsilon() < std::numeric_limits<double>::epsilon();

template <typename scalar>
struct containers {
  using R = std::vector<scalar>;
  using Q = std::deque<scalar>;
  using P = std::map<size_t, scalar>;
};

//! A symmetric matrix with well separated low-lying eigenvalues
template <typename scalar>
struct DiagonallyDominantProblem : molpro::linalg::itsolv::Problem<typename containers<scalar>::R> {
  using R = typename containers<scalar>::R;

  size_t n;
  explicit DiagonallyDominantProblem(size_t n) : n(n) {}

  scalar matrix(size_t i, size_t j) const { return i == j ? scalar(i + 1) : scalar(0.001L) * scalar(i + j); }

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
  std::vector<long double> eigenvalues;
  long double largest_error;
};

template <typename scalar>
Outcome solve(size_t n, size_t nroot, long double convergence_threshold) {
  using C = containers<scalar>;
  auto solver =
      molpro::linalg::itsolv::create_LinearEigensystem<typename C::R, typename C::Q, typename C::P>("Davidson");
  solver->set_n_roots(nroot);
  solver->set_convergence_threshold(scalar(convergence_threshold));
  solver->set_verbosity(0);
  DiagonallyDominantProblem<scalar> problem(n);
  std::vector<typename C::R> params(nroot, typename C::R(n, 0)), actions(nroot, typename C::R(n, 0));
  for (size_t i = 0; i < nroot; ++i)
    params[i][i] = 1; // initial guess: unit vectors on the lowest diagonal elements
  Outcome outcome{};
  outcome.converged = solver->solve(params, actions, problem);
  for (size_t i = 0; i < nroot; ++i)
    outcome.eigenvalues.push_back(static_cast<long double>(solver->eigenvalues()[i]));
  outcome.largest_error = 0;
  for (auto e : solver->errors())
    outcome.largest_error = std::max(outcome.largest_error, static_cast<long double>(e));
  return outcome;
}

constexpr size_t n = 20;
constexpr size_t nroot = 3;

} // namespace

template class molpro::linalg::itsolv::SolverFactory<std::vector<long double>, std::deque<long double>,
                                                     std::map<size_t, long double>>;

//! The double-precision run is the reference: it must converge at a threshold double can reach
TEST(precision_e2e, davidson_double) {
  const auto outcome = solve<double>(n, nroot, 1e-13L);
  EXPECT_TRUE(outcome.converged);
  EXPECT_LT(outcome.largest_error, 1e-13L);
}

//! The same solver in extended precision converges below the double-precision residual floor
TEST(precision_e2e, davidson_long_double) {
  if constexpr (!long_double_is_extended)
    GTEST_SKIP() << "long double is not wider than double on this platform";
  const auto outcome = solve<long double>(n, nroot, 1e-17L);
  EXPECT_TRUE(outcome.converged);
  EXPECT_LT(outcome.largest_error, 1e-17L);
}

namespace {

/*!
 * @brief Minimise 1/2 (x-1)^T A (x-1), whose solution is the vector of ones.
 *
 * Exercises the non-linear solvers, and with them the BFGS line search -- whose interpolant is the
 * other place the library used to compute in double regardless of the surrounding precision.
 */
template <typename scalar>
struct QuadraticForm : molpro::linalg::itsolv::Problem<typename containers<scalar>::R> {
  using R = typename containers<scalar>::R;

  size_t n;
  explicit QuadraticForm(size_t n) : n(n) {}

  scalar matrix(size_t i, size_t j) const { return i == j ? scalar(i + 1) : scalar(0.001L) * scalar(i + j); }

  bool diagonals(R& d) const override {
    for (size_t i = 0; i < d.size(); i++)
      d[i] = matrix(i, i);
    return true;
  }

  scalar residual(const R& v, R& a) const override {
    scalar value = 0;
    for (size_t i = 0; i < a.size(); i++) {
      a[i] = 0;
      for (size_t j = 0; j < a.size(); j++)
        a[i] += matrix(i, j) * (v[j] - 1);
      value += scalar(0.5) * a[i] * (v[i] - 1);
    }
    return value;
  }
};

//! \returns the largest deviation of the converged solution from the vector of ones
template <typename scalar, typename Factory>
long double solve_nonlinear(size_t n, long double convergence_threshold, Factory factory) {
  using C = containers<scalar>;
  auto solver = factory();
  solver->set_convergence_threshold(scalar(convergence_threshold));
  solver->set_verbosity(0);
  QuadraticForm<scalar> problem(n);
  typename C::R x(n, 0), residual(n, 0);
  EXPECT_TRUE(solver->solve(x, residual, problem));
  long double error = 0;
  for (size_t i = 0; i < n; ++i)
    error = std::max(error, std::abs(static_cast<long double>(x[i]) - 1));
  return error;
}

} // namespace

//! The DIIS solver for non-linear equations, in both precisions
TEST(precision_e2e, diis_nonlinear_equations) {
  if constexpr (!long_double_is_extended)
    GTEST_SKIP() << "long double is not wider than double on this platform";
  using Cd = containers<double>;
  using Cl = containers<long double>;
  const auto in_double = solve_nonlinear<double>(10, 1e-12L, [] {
    return molpro::linalg::itsolv::create_NonLinearEquations<typename Cd::R, typename Cd::Q, typename Cd::P>("DIIS");
  });
  const auto in_long_double = solve_nonlinear<long double>(10, 1e-16L, [] {
    return molpro::linalg::itsolv::create_NonLinearEquations<typename Cl::R, typename Cl::Q, typename Cl::P>("DIIS");
  });
  EXPECT_LT(in_double, 1e-11L);
  // below the double-precision floor: only reachable if the whole solve really ran in extended precision
  EXPECT_LT(in_long_double, 1e-15L);
  EXPECT_LT(in_long_double, in_double);
}

//! The BFGS optimiser, whose line search interpolates in the working precision
TEST(precision_e2e, bfgs_optimize) {
  if constexpr (!long_double_is_extended)
    GTEST_SKIP() << "long double is not wider than double on this platform";
  using Cd = containers<double>;
  using Cl = containers<long double>;
  const auto in_double = solve_nonlinear<double>(10, 1e-12L, [] {
    return molpro::linalg::itsolv::create_Optimize<typename Cd::R, typename Cd::Q, typename Cd::P>("BFGS");
  });
  const auto in_long_double = solve_nonlinear<long double>(10, 1e-16L, [] {
    return molpro::linalg::itsolv::create_Optimize<typename Cl::R, typename Cl::Q, typename Cl::P>("BFGS");
  });
  EXPECT_LT(in_double, 1e-11L);
  EXPECT_LT(in_long_double, 1e-15L);
  EXPECT_LT(in_long_double, in_double);
}

//! Both precisions must agree on the answer, to within what double can resolve
TEST(precision_e2e, davidson_precisions_agree) {
  // this comparison is meaningful at any precision, including where the two are the same type
  const auto in_double = solve<double>(n, nroot, 1e-13L);
  const auto in_long_double = solve<long double>(n, nroot, long_double_is_extended ? 1e-17L : 1e-13L);
  ASSERT_EQ(in_double.eigenvalues.size(), in_long_double.eigenvalues.size());
  for (size_t i = 0; i < in_double.eigenvalues.size(); ++i)
    EXPECT_NEAR(double(in_double.eigenvalues[i]), double(in_long_double.eigenvalues[i]), 1e-11);
}

namespace {

//! A non-symmetric matrix whose second and third eigenvalue are a complex pair
template <typename scalar>
struct ComplexPairProblem : molpro::linalg::itsolv::Problem<typename containers<scalar>::R> {
  using R = typename containers<scalar>::R;

  size_t n;
  explicit ComplexPairProblem(size_t n) : n(n) {}

  scalar matrix(size_t i, size_t j) const {
    if (i == 1 and j == 2)
      return 1;
    if (i == 2 and j == 1)
      return -1;
    return i == j ? scalar(i + 1) : scalar(0.01L) / scalar(1 + (i > j ? i - j : j - i));
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

} // namespace

//! Both members of a complex pair have the same error, the norm of their complex residual, in long double too
TEST(precision_e2e, davidson_long_double_complex_pair_errors) {
  using C = containers<long double>;
  auto solver = molpro::linalg::itsolv::create_LinearEigensystem<C::R, C::Q, C::P>("Davidson");
  solver->set_n_roots(3);
  solver->set_hermiticity(false);
  solver->set_convergence_threshold(1e-10L);
  solver->set_verbosity(0);
  ComplexPairProblem<long double> problem(n);
  std::vector<C::R> params(3, C::R(n, 0)), actions(3, C::R(n, 0));
  ASSERT_TRUE(solver->solve(params, actions, problem, true));
  solver->solution({0, 1, 2}, params, actions);
  long double r_real = 0, r_imag = 0;
  for (size_t i = 0; i < n; i++) {
    r_real += actions[1][i] * actions[1][i];
    r_imag += actions[2][i] * actions[2][i];
  }
  const auto errors = solver->errors();
  EXPECT_EQ(errors.at(1), errors.at(2));
  EXPECT_NEAR(errors.at(1), std::sqrt(r_real + r_imag), 1e-17L);
}
