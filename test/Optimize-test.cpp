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
      hmat(i, i) = (i + 2) * param;
    //    std::cout << "hmat " << hmat << std::endl;
  }

  double action(const Rvector& psx, Rvector& outputs) {
    // f = 0.5 x.h.x/x.x + lambda (x.x-1)^2
    // df/dx = x.x^{-1} (h x - x.h.x x.x^{-1} x) + 4 lambda (x.x-1) x
    // f = 0.5 x.h.x/x.x + 0.5 * lambda (x.x0-1)^2
    // df/dx = x.x^{-1} (h x - x.h.x x.x^{-1} x) + lambda (x.x0-1) x0
    // f = 0.5 * (x-b).h.(x-b) where b=[111...]
    // df/dx = h x - h.b

    const double lambda = 1e0;
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


  auto set_options(std::shared_ptr<IOptimize<Rvector, Qvector, Pvector>>& solver, std::shared_ptr<Logger>& logger) {
    auto options = CastOptions::Optimize(solver->get_options());
    options->convergence_threshold = 1.0e-8;
    //    options->norm_thresh = 1.0e-14;
    //    options->svd_thresh = 1.0e-10;
    options->n_roots = 1;
    options->max_size_qspace = 10;
    solver->set_options(*options);
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

  void test_quadratic_form(const std::string& title = "", const int n_working_vectors_max = 0) {
    int nroot = 1;
    {
      int np = 0;
      {
        molpro::cout << "\n\n*** " << title << ",  problem dimension " << n
                     << ", n_working_vectors_max = " << n_working_vectors_max << std::endl;
        auto [solver, logger] = molpro::test::create_Optimize();
        auto options = set_options(solver, logger);
        int nwork = 1;
        // Create initial subspace. This is iteration 0.
        Rvector x(n, 0), g(n);
        x.front() = 1;
        size_t n_iter = 1;
        for (auto iter = 1; iter < 1000; iter++, ++n_iter) {
          auto value = action(x, g);
          auto precon = solver->add_value(x, value, g);
          if (precon)
            update(g);
          if (verbosity > 0)
            solver->report();
          nwork = solver->end_iteration(x, g);
          if (verbosity > 1)
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
        EXPECT_NEAR(solver->value(), 0, 2e-9);
        const auto nR_creations = solver->statistics().r_creations;
        if (verbosity > 0)
          std::cout << "R creations = " << nR_creations << std::endl;
        EXPECT_LE(nR_creations, (nroot + 1) * n_iter);
        std::vector<double> parameters(n), residual(n);
        std::vector<int> roots(1, 0);
        solver->solution(parameters, residual);
        EXPECT_LE(std::sqrt(dot(residual, residual)), options->convergence_threshold.value());
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
    test_quadratic_form(std::to_string(n) + "/" + std::to_string(param));
  }
}
