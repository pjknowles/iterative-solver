#include "test.h"
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
struct RSPT : ::testing::Test {
  using MatrixXdr = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using MatrixXdc = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
  using vectorP = std::vector<scalar>;

  size_t n = 0;
  size_t verbosity = 0;
  MatrixXdc hmat;
  Eigen::VectorXd h0;

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
    h0.resize(n);
    std::ifstream f0(std::string{"./"} + file + ".h0");
    for (auto i = 0; i < n; i++)
      f0 >> h0(i);
    std::cout << "h0" << h0.adjoint() << std::endl;
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

  void update(Rvector& psg) {
    double e0 = 1e50;
    for (size_t i = 0; i < n; i++)
      if (h0(i) < e0)
        e0 = h0(i);
    for (size_t i = 0; i < n; i++) {
      psg[i] = -psg[i] / (1e-12 - e0 + h0(i));
    }
  }

  auto initial_guess(Rvector& x) {
    std::vector<double> diagonals;
    diagonals.reserve(n);
    for (auto i = 0; i < n; i++)
      diagonals.push_back(h0(i));
    x[(std::min_element(diagonals.begin(), diagonals.end()) - diagonals.begin())] = 1; // initial guess
  }

  void test_eigen(const std::string& title = "", const int n_working_vectors_max = 0) {
    {
      molpro::cout << "\n\n*** " << title << " by perturbation theory" << std::endl;
      auto solver = molpro::linalg::itsolv::create_LinearEigensystem<Rvector, Qvector, Pvector>("RSPT");
      int nwork = 1;
      Rvector x(n, 0);
      Rvector g(n, 0);
      initial_guess(x);
      for (auto iter = 1; iter < 5; iter++) {
        action(x, g);
        std::cout << "x " << x << std::endl;
        std::cout << "g " << g << std::endl;
        nwork = solver->add_vector(x, g);
        std::cout << "g after add_vector " << g << std::endl;
        update(g); // TODO e0
        std::cout << "g after update " << g << std::endl;
        solver->report();
        nwork = solver->end_iteration(x, g);
      }
    }
  }
};

TEST_F(RSPT, file_eigen) {
  for (const auto& file : std::vector<std::string>{"bh", "hf"}) {
    load_matrix(file, file == "phenol" ? 0 : 1e-8);
    test_eigen(file);
  }
}
