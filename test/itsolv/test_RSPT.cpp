#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <molpro/iostream.h>
#include <numeric>
#include <regex>
#include <vector>

#include "vector_types.h"
#include <molpro/linalg/itsolv/CastOptions.h>
#include <molpro/linalg/itsolv/SolverFactory.h>
#include <molpro/linalg/itsolv/helper.h>

// Find lowest eigensolution of a matrix obtained from an external file
// Storage of vectors in-memory via class Rvector

using molpro::linalg::array::Span;
using molpro::linalg::itsolv::CastOptions;
using molpro::linalg::itsolv::CVecRef;
using molpro::linalg::itsolv::cwrap;
using molpro::linalg::itsolv::VecRef;
using molpro::linalg::itsolv::wrap;
struct LinearEigensystemF : ::testing::Test {
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

  void action(const Rvector& psx, Rvector& outputs) {
      auto x = Eigen::Map<const Eigen::VectorXd>(psx.data(), psx.size());
      auto r = Eigen::Map<Eigen::VectorXd>(outputs.data(), psx.size());
      r.fill(0);
      r += hmat * x;
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

  auto initial_guess(Rvector& x) {
    std::vector<double> diagonals;
    diagonals.reserve(n);
    for (auto i = 0; i < n; i++)
      diagonals.push_back(hmat(i, i));
    x[(std::min_element(diagonals.begin(), diagonals.end()) - diagonals.begin())] = 1; // initial guess
  }


  void test_eigen(const std::string& title = "", const int n_working_vectors_max = 0) {
    bool hermitian = check_mat_hermiticity();
    auto [expected_eigensolutions, expected_eigenvalues] = solve_full_problem(hermitian);
    {
      molpro::cout << "\n\n*** " << title << " by perturbation theory" << std::endl;
      auto solver = molpro::linalg::itsolv::create_LinearEigensystem<Rvector, Qvector, Pvector>("RSPT");
      int nwork = 1;
      Rvector x(n,0);
      Rvector g(n,0);
      initial_guess(x);
      for (auto iter = 1; iter < 5; iter++) {
        action(x, g);
        std::cout << "x "<<x<<std::endl;
        std::cout << "g "<<g<<std::endl;
        nwork = solver->add_vector(x, g);
        std::cout << "g after add_vector "<<g<<std::endl;
        update(g, solver->working_set_eigenvalues());
        solver->report();
        nwork = solver->end_iteration(x, g);
      }
    }
  }
};

TEST_F(LinearEigensystemF, file_eigen) {
  for (const auto& file : std::vector<std::string>{"bh", "hf", "phenol"}) {
    load_matrix(file, file == "phenol" ? 0 : 1e-8);
    test_eigen(file);
  }
}

TEST_F(LinearEigensystemF, n_eigen) {
  size_t n = 100;
  double param = 1;
  //  for (auto param : std::vector<double>{.01, .1, 1, 10, 100}) {
  //  for (auto param : std::vector<double>{.01, .1, 1}) {
  for (auto param : std::vector<double>{1}) {
    load_matrix(n, "", param);
    test_eigen(std::to_string(n) + "/" + std::to_string(param));
  }
}
TEST_F(LinearEigensystemF, nonhermitian_eigen) {
  size_t n = 30;
  double param = 1;
  //  for (auto param : std::vector<double>{.01, .1, 1, 10, 100}) {
  //  for (auto param : std::vector<double>{.01, .1, 1}) {
  for (auto param : std::vector<double>{1}) {
    load_matrix(n, "", param, false);
    test_eigen(std::to_string(n) + "/" + std::to_string(param) + ", non-hermitian");
  }
}

TEST_F(LinearEigensystemF, small_eigen) {
  for (int n = 1; n < 20; n++) {
    double param = 1;
    load_matrix(n, "", param);
    test_eigen(std::to_string(n) + "/" + std::to_string(param));
    test_eigen(std::to_string(n) + "/" + std::to_string(param), 1);
  }
}

TEST_F(LinearEigensystemF, symmetry_eigen) {
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
