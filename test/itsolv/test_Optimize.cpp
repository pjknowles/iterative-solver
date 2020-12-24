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
#include <molpro/linalg/itsolv/Optimize.h>
#include <molpro/linalg/itsolv/helper.h>
#include <molpro/linalg/itsolv/subspace/SubspaceSolverOptBFGS.h>
#include <molpro/linalg/itsolv/subspace/SubspaceSolverOptSD.h>
#include <molpro/linalg/itsolv/SolverFactory.h>
#include "vector_types.h"
using molpro::linalg::itsolv::CastOptions;
using molpro::linalg::itsolv::Logger;
using molpro::linalg::array::Span;
using molpro::linalg::itsolv::CVecRef;
using molpro::linalg::itsolv::cwrap;
using molpro::linalg::itsolv::VecRef;
using molpro::linalg::itsolv::wrap;
struct OptimizeF : ::testing::Test {
  using scalar = double;
  using MatrixXdc = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

  size_t n = 0;
  size_t verbosity = 0;
  MatrixXdc hmat;

  void load_matrix(int dimension, const std::string& type = "", double param = 1, bool hermitian = true) {
    n = dimension;
    hmat.resize(n, n);
    hmat.fill(1);
    for (int i = 0; i < n; i++)
      hmat(i, i) = (i + 2) * param;
    //    std::cout << "hmat " << hmat << std::endl;
  }

  double action(const Rvector& psx, Rvector& outputs) {
    // f = 0.5 * (x-b).h.(x-b) where b=[111...]
    // df/dx = h x - h.b
    auto x = Eigen::Map<const Eigen::VectorXd>(psx.data(), psx.size());
    auto r = Eigen::Map<Eigen::VectorXd>(outputs.data(), psx.size());
    Eigen::VectorXd b;
    b.resize(x.size());
    b.fill(1);
    r = hmat * (x - b);
    double value = 0.5 * (x - b).dot(hmat * (x - b));
    //    std::cout << "action, x:\n" << x.transpose() << std::endl;
    //    std::cout << "action, r:\n" << r.transpose() << std::endl;
    //    std::cout << "action, value: " << value << std::endl;
    //    std::cout << "action, norm squared: " << x.dot(x) << std::endl;
    return value;
  }
  template <class scalar>
  scalar dot(const std::vector<scalar>& a, const std::vector<scalar>& b) {
    return std::inner_product(std::begin(a), std::end(a), std::begin(b), scalar(0));
  }

  void update(Rvector& psg) {
    for (size_t i = 0; i < n; i++)
      psg[i] = -psg[i] / hmat(i, i);
  }

  void test_quadratic_form(const std::string& method, const std::string& title = "", const int n_working_vectors_max = 0) {
    int nroot = 1;
    {
      int np = 0;
      {
        molpro::cout << "\n\n*** " << title << ",  problem dimension " << n
                     << ", n_working_vectors_max = " << n_working_vectors_max << std::endl;

        auto solver = molpro::linalg::itsolv::create_Optimize<Rvector ,Qvector >(method,"convergence_threshold=1e-8,n_roots=1,max_size_qspace=10");
        int nwork = 1;
        // Create initial subspace. This is iteration 0.
        Rvector x(n, 0), g(n);
        x.front() = 1;
        size_t n_iter = 1;
        for (auto iter = 1; iter < 1000 && nwork > 0; iter++, ++n_iter) {
          auto value = action(x, g);
          if (solver->add_value(x, value, g))
            update(g);
          if (verbosity > 0)
            solver->report();
          nwork = solver->end_iteration(x, g);
        }
        std::cout << "Error={ ";
        for (const auto& e : solver->errors())
          std::cout << e << " ";
        std::cout << "} after " << solver->statistics().iterations << " iterations, "
                  << solver->statistics().r_creations << " R vectors" << std::endl;
        EXPECT_THAT(solver->errors(), ::testing::Pointwise(::testing::DoubleNear(2 * solver->convergence_threshold()),
                                                           std::vector<double>(nroot, double(0))));
        EXPECT_NEAR(solver->value(), 0, 2e-9);
        const auto nR_creations = solver->statistics().r_creations;
        if (verbosity > 0)
          std::cout << "R creations = " << nR_creations << std::endl;
        EXPECT_LE(nR_creations, (nroot + 1) * n_iter);
        std::vector<double> parameters(n), residual(n);
        std::vector<int> roots(1, 0);
        solver->solution(parameters, residual);
        EXPECT_LE(std::sqrt(dot(residual, residual)), solver->convergence_threshold());
        EXPECT_THAT(parameters, ::testing::Pointwise(::testing::DoubleNear(solver->convergence_threshold()),
                                                     std::vector<double>(n, double(1))));
      }
    }
  }
};

TEST_F(OptimizeF, small_quadratic_form) {
  for (int n = 1; n < 51; n++) {
    double param = 10;
    load_matrix(n, "", param);
    test_quadratic_form("BFGS",std::to_string(n) + "/" + std::to_string(param));
    test_quadratic_form("SD", std::to_string(n) + "/" + std::to_string(param));
  }
}
