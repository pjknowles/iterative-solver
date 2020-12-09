#include "create_solver.h"
#include "test.h"

#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <molpro/iostream.h>
#include <numeric>
#include <regex>
#include <vector>

#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/helper.h>

// Find lowest eigensolution of a matrix obtained from an external file
// Storage of vectors in-memory via class Rvector

using molpro::linalg::array::Span;
using molpro::linalg::itsolv::CVecRef;
using molpro::linalg::itsolv::cwrap;
using molpro::linalg::itsolv::VecRef;
using molpro::linalg::itsolv::wrap;
struct OptimizeF : ::testing::Test {
  using scalar = double;
  using MatrixXdr = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using MatrixXdc = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
  using vectorP = std::vector<scalar>;

  size_t n = 0;
  size_t verbosity = 0;
  MatrixXdc hmat;

  void load_matrix(int dimension, const std::string& type = "", double param = 1, bool hermitian = true) {
    n = dimension;
    hmat.resize(n, n);
    hmat.fill(1);
    for (int i = 0; i < n; i++)
      hmat(i, i) = i * param;
    if (not hermitian)
      for (int i = 0; i < n; i++)
        for (int j = 0; j < i; j++)
          hmat(i, j) *= .95;
  }
  void load_matrix(const std::string& file, double degeneracy_split = 0) {
    std::ifstream f(std::string{"./"} + file + ".hamiltonian");
    f >> n;
    hmat.resize(n, n);
    for (auto i = 0; i < n; i++)
      for (auto j = 0; j < n; j++)
        f >> hmat(i, j);
    // split degeneracies
    for (auto i = 0; i < n; i++)
      hmat(i, i) += degeneracy_split * i;
  }

  double action(const std::vector<Rvector>& psx, std::vector<Rvector>& outputs) {
    size_t k = 0;
    auto x = Eigen::Map<const Eigen::VectorXd>(psx.at(k).data(), psx.at(k).size());
    auto r = Eigen::Map<Eigen::VectorXd>(outputs.at(k).data(), psx.at(k).size());
    r.fill(0);
    r += 2 * hmat * x;
    double value = 0.5 * x.dot(r) / x.dot(x);
    r -= 2 * value * x;
    return value;
  }
  template <class scalar>
  scalar dot(const std::vector<scalar>& a, const std::vector<scalar>& b) {
    return std::inner_product(std::begin(a), std::end(a), std::begin(b), scalar(0));
  }
  void residual(const std::vector<Rvector>& psx, std::vector<Rvector>& actions, const std::vector<double>& evals) {
    action(psx, actions);
    for (size_t k = 0; k < psx.size(); k++) {
      auto x = Eigen::Map<const Eigen::VectorXd>(psx.at(k).data(), psx.at(k).size());
      auto r = Eigen::Map<Eigen::VectorXd>(actions.at(k).data(), psx.at(k).size());
      r -= evals[k] * x;
    }
  }

  void update(std::vector<Rvector>& psg, std::vector<scalar> shift = std::vector<scalar>()) {
    for (size_t k = 0; k < psg.size() && k < shift.size(); k++)
      for (size_t i = 0; i < n; i++) {
        psg[k][i] = -psg[k][i] / (1e-12 - shift[k] + hmat(i, i));
      }
  }

  std::tuple<std::map<double, std::vector<double>>, std::vector<double>> solve_full_problem(bool hermitian) {
    std::map<double, std::vector<double>> expected_eigensolutions;
    std::vector<double> expected_eigenvalues, eigenvector;
    std::vector<double> hmat_row(n * n);
    std::vector<double> metric(n * n);
    for (size_t i = 0; i < n; ++i)
      metric[i * (n + 1)] = 1;
    for (size_t i = 0, ij = 0; i < n; ++i)
      for (size_t j = 0; j < n; ++j, ++ij)
        hmat_row[ij] = hmat(i, j);
    molpro::linalg::itsolv::eigenproblem(eigenvector, expected_eigenvalues, hmat_row, metric, n, hermitian, 1.0e-14, 0);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        expected_eigensolutions[expected_eigenvalues[i]].push_back(
            eigenvector[n * i + j]); // won't work for degenerate eigenvalues!
    std::sort(expected_eigenvalues.begin(), expected_eigenvalues.end());
    return std::make_tuple(expected_eigensolutions, expected_eigenvalues);
  }

  // testing the test: that sort of eigenvalues is the same thing as std::map sorting of keys
  void check_eigenvectors_map(const std::map<double, std::vector<double>>& expected_eigensolutions,
                              const std::vector<double>& expected_eigenvalues) {
    std::vector<double> testing;
    for (const auto& ee : expected_eigensolutions)
      testing.push_back(ee.first);
    ASSERT_THAT(expected_eigenvalues, ::testing::Pointwise(::testing::DoubleEq(), testing));
  }

  auto create_containers(const size_t nW) {
    std::vector<Rvector> g;
    std::vector<Rvector> x;
    for (size_t i = 0; i < nW; i++) {
      x.emplace_back(n, 0);
      g.emplace_back(n);
    }
    return std::make_tuple(x, g);
  }

  auto initial_guess(std::vector<Rvector>& x) {
    const auto n_roots = x.size();
    std::vector<size_t> guess;
    std::vector<double> diagonals;
    for (auto i = 0; i < n; i++) {
      diagonals.push_back(hmat(i, i));
    }
    for (size_t root = 0; root < n_roots; root++) {
      guess.push_back(std::min_element(diagonals.begin(), diagonals.end()) - diagonals.begin()); // initial guess
      *std::min_element(diagonals.begin(), diagonals.end()) = 1e99;
      x[root][guess.back()] = 1; // initial guess
    }
  }

  void apply_p(const std::vector<vectorP>& pvectors, const CVecRef<Pvector>& pspace, const VecRef<Rvector>& action) {
    for (size_t i = 0; i < pvectors.size(); i++) {
      auto& actioni = (action[i]).get();
      for (size_t pi = 0; pi < pspace.size(); pi++) {
        const auto& p = pspace[pi].get();
        for (const auto& pel : p)
          for (size_t j = 0; j < n; j++)
            actioni[j] += hmat(j, pel.first) * pel.second * pvectors[i][pi];
      }
    }
  }

  auto initial_pspace(const size_t np) {
    std::vector<scalar> PP;
    std::vector<Pvector> pspace;
    if (np > 0) {
      std::vector<double> diagonals;
      for (auto i = 0; i < n; i++)
        diagonals.push_back(hmat(i, i));
      for (size_t p = 0; p < np; p++) {
        std::map<size_t, scalar> pp;
        pp.insert({std::min_element(diagonals.begin(), diagonals.end()) - diagonals.begin(), 1});
        pspace.push_back(pp);
        *std::min_element(diagonals.begin(), diagonals.end()) = 1e99;
      }
      PP.reserve(np * np);
      for (const auto& i : pspace)
        for (const auto& j : pspace) {
          const size_t kI = i.begin()->first;
          const size_t kJ = j.begin()->first;
          PP.push_back(hmat(kI, kJ));
        }
    }
    return std::make_tuple(pspace, PP);
  }

  auto set_options(std::shared_ptr<IOptimize<Rvector, Qvector, Pvector>>& solver, std::shared_ptr<Logger>& logger) {
    auto options = CastOptions::Optimize(solver->get_options());
    options->convergence_threshold = 1.0e-8;
    //    options->norm_thresh = 1.0e-14;
    //    options->svd_thresh = 1.0e-10;
    options->max_size_qspace = 10;
    solver->set_options(options);
    options = CastOptions::Optimize(solver->get_options());
    molpro::cout << "convergence threshold = " << options->convergence_threshold.value()
                 << ", svd thresh = " << options->svd_thresh.value()
                 << ", norm thresh = " << options->norm_thresh.value()
                 << ", max size of Q = " << options->max_size_qspace.value() << std::endl;
    logger->max_trace_level = molpro::linalg::itsolv::Logger::None;
    logger->max_warn_level = molpro::linalg::itsolv::Logger::Error;
    logger->data_dump = false;
    return options;
  }

  void test_eigen(const std::string& title = "", const int n_working_vectors_max = 0) {
    auto d = (hmat - hmat.transpose()).norm();
    bool hermitian = d < 1e-10;
    auto [expected_eigensolutions, expected_eigenvalues] = solve_full_problem(hermitian);
    check_eigenvectors_map(expected_eigensolutions, expected_eigenvalues);
    if (verbosity > 0)
      std::cout << "expected eigenvalues " << expected_eigenvalues << std::endl;
    int nroot = 1;
    {
      int np = 0;
      {
        molpro::cout << "\n\n*** " << title << ", " << nroot << " roots, problem dimension " << n
                     << ", pspace dimension " << np << ", n_working_vectors_max = " << n_working_vectors_max
                     << std::endl;
        auto [solver, logger] = molpro::test::create_Optimize();
        auto options = set_options(solver, logger);
        int nwork = 1;
        // Create initial subspace. This is iteration 0.
        auto [x, g] = create_containers(1);
        initial_guess(x);
        size_t n_iter = 1;
        for (auto iter = 1; iter < 100; iter++, ++n_iter) {
          auto value = action(x, g);
          auto precon = solver->add_value(x.front(), value, g.front());
          if (precon)
            update(g, {value});
          if (verbosity > 0)
            solver->report();
          nwork = solver->end_iteration(x, g);
          if (verbosity > 0)
            std::cout << "solver.end_iteration returns nwork=" << nwork << std::endl;
          if (nwork == 0)
            break;
        }
        if (verbosity > 0)
          solver->report();
        std::cout << "Error={ ";
        for (const auto& e : solver->errors())
          std::cout << e << " ";
        std::cout << "} after " << solver->statistics().iterations << " iterations, "
                  << solver->statistics().r_creations << " R vectors" << std::endl;
        EXPECT_THAT(solver->errors(), ::testing::Pointwise(::testing::DoubleNear(2 * solver->convergence_threshold()),
                                                           std::vector<double>(nroot, double(0))));
        if (verbosity > 0) {
          std::cout << "expected eigenvalues "
                    << std::vector<double>(expected_eigenvalues.data(), expected_eigenvalues.data() + solver->n_roots())
                    << std::endl;
          std::cout << "obtained eigenvalues " << solver->value() << std::endl;
        }
        EXPECT_NEAR(solver->value(), expected_eigenvalues.front(), 2e-9);
        const auto nR_creations = solver->statistics().r_creations;
        if (verbosity > 0)
          std::cout << "R creations = " << nR_creations << std::endl;
        EXPECT_LE(nR_creations, (nroot + 1) * n_iter);
        std::vector<std::vector<double>> parameters, residuals;
        std::vector<int> roots;
        for (int root = 0; root < solver->n_roots(); root++) {
          parameters.emplace_back(n);
          residuals.emplace_back(n);
          roots.push_back(root);
        }
        solver->solution(roots, parameters, residuals);
        residual(parameters, residuals, {solver->value()});
        for (const auto& r : residuals)
          EXPECT_LE(std::sqrt(dot(r, r)), options->convergence_threshold.value());
        int root = 0;
        for (const auto& ee : expected_eigensolutions) {
          EXPECT_NEAR(ee.first, solver->value(), 1e-10);
          auto overlap_with_reference = std::inner_product(parameters.at(root).begin(), parameters.at(root).end(),
                                                           ee.second.begin(), 0., std::plus(), std::multiplies());
          EXPECT_NEAR(std::abs(overlap_with_reference), 1., options->convergence_threshold.value());
          root++;
          if (root >= solver->n_roots())
            break;
        }
      }
    }
  }
};

TEST_F(OptimizeF, file_eigen) {
  for (const auto& file : std::vector<std::string>{"phenol", "bh", "hf"}) {
    load_matrix(file, file == "phenol" ? 0 : 1e-8);
    test_eigen(file);
  }
}

TEST_F(OptimizeF, n_eigen) {
  size_t n = 100;
  double param = 1;
  //  for (auto param : std::vector<double>{.01, .1, 1, 10, 100}) {
  //  for (auto param : std::vector<double>{.01, .1, 1}) {
  for (auto param : std::vector<double>{1}) {
    load_matrix(n, "", param);
    test_eigen(std::to_string(n) + "/" + std::to_string(param));
  }
}
TEST_F(OptimizeF, nonhermitian_eigen) {
  size_t n = 30;
  double param = 1;
  //  for (auto param : std::vector<double>{.01, .1, 1, 10, 100}) {
  //  for (auto param : std::vector<double>{.01, .1, 1}) {
  for (auto param : std::vector<double>{1}) {
    load_matrix(n, "", param, false);
    test_eigen(std::to_string(n) + "/" + std::to_string(param) + ", non-hermitian");
  }
}

TEST_F(OptimizeF, small_eigen) {
  for (int n = 1; n < 20; n++) {
    double param = 1;
    load_matrix(n, "", param);
    test_eigen(std::to_string(n) + "/" + std::to_string(param));
    test_eigen(std::to_string(n) + "/" + std::to_string(param), 1);
  }
}

TEST_F(OptimizeF, symmetry_eigen) {
  for (int n = 1; n < 10; n++) {
    double param = 1;
    load_matrix(n, "", param);
    for (auto i = 0; i < n; i++)
      for (auto j = 0; j < n; j++)
        if (((i % 3) == 0 and (j % 3) != 0) or ((j % 3) == 0 and (i % 3) != 0))
          hmat(j, i) = 0;
    if (false) {
      std::cout << "hmat " << std::endl;
      for (auto i = 0; i < n; i++) {
        std::cout << "i%3" << i % 3 << std::endl;
        for (auto j = 0; j < n; j++)
          std::cout << " " << hmat(j, i);
        std::cout << std::endl;
      }
    }
    test_eigen(std::to_string(n) + "/sym/" + std::to_string(param));
  }
}
