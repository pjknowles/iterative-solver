#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/array/DistrArrayFile.h>
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/SolverFactory.h>
#include <molpro/mpi.h>

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <complex>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace {
//! MPI is initialised by the tests that need it, and finalised once, after all tests have run
class MPIEnvironment : public ::testing::Environment {
public:
  void TearDown() override {
#ifdef HAVE_MPI_H
    int initialized, finalized;
    MPI_Initialized(&initialized);
    MPI_Finalized(&finalized);
    if (initialized and not finalized)
      molpro::mpi::finalize();
#endif
  }
};
[[maybe_unused]] auto *const mpi_environment = ::testing::AddGlobalTestEnvironment(new MPIEnvironment);
} // namespace

template <typename Container = std::vector<double>>
class ComplexRootsProblem : public molpro::linalg::itsolv::Problem<Container> {
public:
  using container_t = Container;
  using value_t = typename Container::value_type;

  // clang-format off
  Eigen::Matrix<value_t, 4, 4, Eigen::RowMajor> matrix{
      { 2, -1, -2, 0},
      {-1,  1, -3, 1},
      { 0,  2,  3,  0},
	  { 3,  2,  1,  0},
  };
  // clang-format on
  using VecT = Eigen::Vector4d;

  std::size_t dim() const { return matrix.cols(); }

  bool diagonals(container_t &d) const override {
    for (std::size_t i = 0; i < dim(); ++i) {
      d.at(i) = matrix(i, i);
    }

    return true;
  }

  void action(const molpro::linalg::itsolv::CVecRef<container_t> &parameters,
              const molpro::linalg::itsolv::VecRef<container_t> &actions) const override {
    for (std::size_t i = 0; i < parameters.size(); i++) {
      Eigen::Map<const VecT> input(parameters.at(i).get().data(), parameters.at(i).get().size());
      Eigen::Map<VecT> output(actions.at(i).get().data(), actions.at(i).get().size());

      output = matrix * input;
    }
  }
};

TEST(LinearEigensystem, Davidson_complex_roots) {
  molpro::mpi::init();

  {
    ComplexRootsProblem<> problem;

    using Rvector = decltype(problem)::container_t;

    auto solver = molpro::linalg::itsolv::create_LinearEigensystem<Rvector>("Davidson");
    solver->set_n_roots(4);
    solver->set_max_iter(100);

    std::vector<Rvector> params;
    std::vector<Rvector> actions;
    for (std::size_t i = 0; i < solver->n_roots(); ++i) {
      params.emplace_back(problem.dim());
      actions.emplace_back(problem.dim());
    }

    if (!solver->solve(params, actions, problem, true)) {
      FAIL() << "Solver did not converge";
    }

    std::vector<int> roots(solver->n_roots());
    std::iota(roots.begin(), roots.end(), 0);

    solver->solution(roots, params, actions);
    auto eigvals = solver->eigenvalues();

    ASSERT_EQ(eigvals.size(), solver->n_roots());
    const double tol = 1e-9;
    ASSERT_THAT(eigvals.at(0), ::testing::DoubleNear(-0.360430299, tol));
    ASSERT_THAT(eigvals.at(1), ::testing::DoubleNear(1.5975220393, tol));
    ASSERT_THAT(eigvals.at(2), ::testing::DoubleNear(2.3814541296, tol));
    ASSERT_THAT(eigvals.at(3), ::testing::DoubleNear(2.3814541296, tol));
  }
}

//! A non-symmetric, diagonally dominant matrix whose third eigenvalue (by real part) is a member of a complex pair
class SplitComplexPairProblem : public molpro::linalg::itsolv::Problem<std::vector<double>> {
public:
  static constexpr size_t n = 30;
  double element(size_t i, size_t j) const {
    if (i == 2 and j == 3)
      return 1.;
    if (i == 3 and j == 2)
      return -1.;
    return i == j ? 1. + i : 0.01 / (1. + std::abs(double(i) - double(j)));
  }

  bool diagonals(std::vector<double> &d) const override {
    for (size_t i = 0; i < n; ++i)
      d[i] = element(i, i);
    return true;
  }

  void action(const molpro::linalg::itsolv::CVecRef<std::vector<double>> &parameters,
              const molpro::linalg::itsolv::VecRef<std::vector<double>> &actions) const override {
    for (size_t k = 0; k < parameters.size(); ++k)
      for (size_t i = 0; i < n; ++i) {
        actions[k].get()[i] = 0;
        for (size_t j = 0; j < n; ++j)
          actions[k].get()[i] += element(i, j) * parameters[k].get()[j];
      }
  }
};

TEST(LinearEigensystem, Davidson_complex_pair_split_by_n_roots) {
  molpro::mpi::init();
  SplitComplexPairProblem problem;
  Eigen::MatrixXd matrix(problem.n, problem.n);
  for (size_t i = 0; i < problem.n; ++i)
    for (size_t j = 0; j < problem.n; ++j)
      matrix(i, j) = problem.element(i, j);
  const Eigen::VectorXcd reference_complex = Eigen::EigenSolver<Eigen::MatrixXd>(matrix).eigenvalues();
  std::vector<std::complex<double>> reference(reference_complex.begin(), reference_complex.end());
  std::ranges::sort(reference, {}, [](const auto &e) { return e.real(); });
  ASSERT_NEAR(std::abs(reference.at(2) - std::conj(reference.at(3))), 0., 1e-12);
  ASSERT_GT(std::abs(reference.at(2).imag()), 0.5);

  auto solver = molpro::linalg::itsolv::create_LinearEigensystem<std::vector<double>>("Davidson");
  solver->set_n_roots(3);
  solver->set_max_iter(100);
  solver->set_hermiticity(false);
  solver->set_convergence_threshold(1e-8);
  auto params = std::vector<std::vector<double>>(solver->n_roots(), std::vector<double>(problem.n));
  auto actions = params;
  ASSERT_TRUE(solver->solve(params, actions, problem, true));
  const auto eigvals = solver->eigenvalues();
  ASSERT_EQ(eigvals.size(), solver->n_roots());
  for (size_t i = 0; i < eigvals.size(); ++i)
    EXPECT_NEAR(eigvals[i], reference[i].real(), 1e-8);
}

//! Records the warnings issued through it
class WarningRecorder : public molpro::linalg::itsolv::Logger {
public:
  mutable std::vector<std::string> warnings;

protected:
  void default_message_handler(std::string_view ctx, molpro::linalg::itsolv::log::Severity severity,
                               molpro::linalg::itsolv::log::Verbosity verbosity, std::string_view msg) const override {
    if (severity == molpro::linalg::itsolv::log::Severity::Warning)
      warnings.emplace_back(msg);
  }
};

//! The user is told that only the real part of the eigenvector of a root with a truncated complex partner is available
TEST(LinearEigensystem, Davidson_complex_pair_split_by_n_roots_warning) {
  molpro::mpi::init();
  SplitComplexPairProblem problem;
  auto solver = molpro::linalg::itsolv::create_LinearEigensystem<std::vector<double>>("Davidson");
  auto logger = std::make_shared<WarningRecorder>();
  solver->set_logger(logger);
  solver->set_n_roots(3);
  solver->set_max_iter(100);
  solver->set_hermiticity(false);
  solver->set_convergence_threshold(1e-8);
  auto params = std::vector<std::vector<double>>(solver->n_roots(), std::vector<double>(problem.n));
  auto actions = params;
  ASSERT_TRUE(solver->solve(params, actions, problem, true));
  EXPECT_THAT(logger->warnings, ::testing::Contains(::testing::HasSubstr(
                                    "partner of root 3 is not among the requested roots, so only the real part")));
}

//! A complex pair split by the number of roots with a reduced subspace, which has to keep the truncated partner
class Davidson_complex_pair_split_reduced_subspace : public ::testing::TestWithParam<std::string> {};

TEST_P(Davidson_complex_pair_split_reduced_subspace, converges) {
  molpro::mpi::init();
  SplitComplexPairProblem problem;
  Eigen::MatrixXd matrix(problem.n, problem.n);
  for (size_t i = 0; i < problem.n; ++i)
    for (size_t j = 0; j < problem.n; ++j)
      matrix(i, j) = problem.element(i, j);
  const Eigen::VectorXcd reference_complex = Eigen::EigenSolver<Eigen::MatrixXd>(matrix).eigenvalues();
  std::vector<std::complex<double>> reference(reference_complex.begin(), reference_complex.end());
  std::ranges::sort(reference, {}, [](const auto &e) { return e.real(); });

  auto solver = molpro::linalg::itsolv::create_LinearEigensystem<std::vector<double>>("Davidson", GetParam());
  solver->set_n_roots(3);
  solver->set_max_iter(100);
  solver->set_hermiticity(false);
  solver->set_convergence_threshold(1e-8);
  auto params = std::vector<std::vector<double>>(solver->n_roots(), std::vector<double>(problem.n));
  auto actions = params;
  ASSERT_TRUE(solver->solve(params, actions, problem, true));
  const auto eigvals = solver->eigenvalues();
  ASSERT_EQ(eigvals.size(), solver->n_roots());
  for (size_t i = 0; i < eigvals.size(); ++i)
    EXPECT_NEAR(eigvals[i], reference[i].real(), 1e-8);
}

INSTANTIATE_TEST_SUITE_P(LinearEigensystem, Davidson_complex_pair_split_reduced_subspace,
                         ::testing::Values("max_size_qspace=4", "contrib_thresh=1e-3", "max_size_qspace=4,reset_D=3"));

//! The two members of a complex pair end up in different batches of the solution when there are fewer parameter
//! vectors than roots
TEST(LinearEigensystem, Davidson_complex_pair_split_by_batches) {
  molpro::mpi::init();
  SplitComplexPairProblem problem;
  Eigen::MatrixXd matrix(problem.n, problem.n);
  for (size_t i = 0; i < problem.n; ++i)
    for (size_t j = 0; j < problem.n; ++j)
      matrix(i, j) = problem.element(i, j);
  const Eigen::VectorXcd reference_complex = Eigen::EigenSolver<Eigen::MatrixXd>(matrix).eigenvalues();
  std::vector<std::complex<double>> reference(reference_complex.begin(), reference_complex.end());
  std::ranges::sort(reference, {}, [](const auto &e) { return e.real(); });
  ASSERT_NEAR(std::abs(reference.at(2) - std::conj(reference.at(3))), 0., 1e-12);

  auto solver = molpro::linalg::itsolv::create_LinearEigensystem<std::vector<double>>("Davidson");
  solver->set_n_roots(4);
  solver->set_max_iter(100);
  solver->set_hermiticity(false);
  solver->set_convergence_threshold(1e-8);
  auto params = std::vector<std::vector<double>>(3, std::vector<double>(problem.n));
  auto actions = params;
  ASSERT_TRUE(solver->solve(params, actions, problem, true));
  const auto eigvals = solver->eigenvalues();
  ASSERT_EQ(eigvals.size(), solver->n_roots());
  for (size_t i = 0; i < eigvals.size(); ++i)
    EXPECT_NEAR(eigvals[i], reference[i].real(), 1e-8);
}

//! A symmetric, diagonally dominant matrix
class SymmetricProblem : public molpro::linalg::itsolv::Problem<std::vector<double>> {
public:
  static constexpr size_t n = 60;
  double element(size_t i, size_t j) const { return i == j ? 1. + i : 0.1 / (1. + std::abs(double(i) - double(j))); }

  bool diagonals(std::vector<double> &d) const override {
    for (size_t i = 0; i < n; ++i)
      d[i] = element(i, i);
    return true;
  }

  void action(const molpro::linalg::itsolv::CVecRef<std::vector<double>> &parameters,
              const molpro::linalg::itsolv::VecRef<std::vector<double>> &actions) const override {
    for (size_t k = 0; k < parameters.size(); ++k)
      for (size_t i = 0; i < n; ++i) {
        actions[k].get()[i] = 0;
        for (size_t j = 0; j < n; ++j)
          actions[k].get()[i] += element(i, j) * parameters[k].get()[j];
      }
  }
};

template <class Qvector>
std::vector<double> solve_symmetric_problem() {
  SymmetricProblem problem;
  auto solver = molpro::linalg::itsolv::create_LinearEigensystem<std::vector<double>, Qvector>("Davidson");
  solver->set_n_roots(3);
  solver->set_max_iter(100);
  solver->set_hermiticity(true);
  solver->set_convergence_threshold(1e-10);
  auto params = std::vector<std::vector<double>>(solver->n_roots(), std::vector<double>(problem.n));
  auto actions = params;
  EXPECT_TRUE(solver->solve(params, actions, problem, true));
  return solver->eigenvalues();
}

// The R vectors are contiguous in memory and the Q vectors on disk, as handled by ArrayHandlerIterableDDisk and
// ArrayHandlerDDiskIterable
TEST(LinearEigensystem, Davidson_vector_R_disk_Q) {
  molpro::mpi::init();
  const auto reference = solve_symmetric_problem<std::vector<double>>();
  const auto result = solve_symmetric_problem<molpro::linalg::array::DistrArrayFile>();
  ASSERT_EQ(result.size(), reference.size());
  for (size_t i = 0; i < reference.size(); ++i)
    EXPECT_NEAR(result[i], reference[i], 1e-9);
}

#include <molpro/linalg/itsolv/SolverFactory-implementation.h>
template class molpro::linalg::itsolv::SolverFactory<std::vector<double>, std::vector<double>,
                                                     std::map<size_t, double>>;
template class molpro::linalg::itsolv::SolverFactory<std::vector<double>, molpro::linalg::array::DistrArrayFile,
                                                     std::map<size_t, double>>;
