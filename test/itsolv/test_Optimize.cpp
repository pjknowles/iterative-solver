#include "test.h"

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <molpro/iostream.h>
#include <numeric>
#include <vector>

#include "vector_types.h"
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/SolverFactory.h>
#ifndef NOFORTRAN
extern "C" int test_optimizef(double* matrix, size_t n);
#endif
struct OptimizeF : ::testing::Test {
  using MatrixXdc = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

  size_t n = 0;
  size_t verbosity = 2;
  MatrixXdc hmat;

  void load_matrix(int dimension, const std::string& type = "", double param = 1, bool hermitian = true) {
    n = dimension;
    hmat.resize(n, n);
    hmat.fill(1);
    for (size_t i = 0; i < n; i++)
      hmat(i, i) = (i + 2) * param;
    //    molpro::cout << "hmat " << hmat << std::endl;
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
    //    molpro::cout << "action, x:\n" << x.transpose() << std::endl;
    //    molpro::cout << "action, r:\n" << r.transpose() << std::endl;
    //    molpro::cout << "action, value: " << value << std::endl;
    //    molpro::cout << "action, norm squared: " << x.dot(x) << std::endl;
    return value;
  }
  template <class scalar>
  scalar dot(const std::vector<scalar>& a, const std::vector<scalar>& b) {
    return std::inner_product(std::begin(a), std::end(a), std::begin(b), scalar(0));
  }

  void update(Rvector& psg) {
    for (size_t i = 0; i < n; i++)
      psg[i] = psg[i] / hmat(i, i);
  }

  void test_quadratic_form(const std::string& method, const std::string& title = "",
                           const int n_working_vectors_max = 0) {
    int nroot = 1;
    {
//      int np = 0;
      {
        molpro::cout << "\n\n*** " << title << ",  problem dimension " << n << ", method = " << method
                     << ", n_working_vectors_max = " << n_working_vectors_max << std::endl;

        auto solver = molpro::linalg::itsolv::create_Optimize<Rvector, Qvector>(
            method, "convergence_threshold=1e-8,max_size_qspace=6");
        int nwork = 1;
        // Create initial subspace. This is iteration 0.
        Rvector x(n, 0), g(n);
        x.front() = 1;
        size_t n_iter = 1;
        for (auto iter = 1; iter < 1000 && nwork > 0; iter++, ++n_iter) {
          auto value = action(x, g);
          //                    molpro::cout << "Iteration "<<iter<<std::endl;
          //                    molpro::cout << "x\n"<<x<< std::endl;
          //                    molpro::cout << "g\n"<<g<< std::endl;
          if (solver->add_vector(x, g, value) > 0) {
            //            molpro::cout << "before update, g\n"<<g<< std::endl;
            update(g);
          }
          //          molpro::cout << "after update, g\n"<<g<< std::endl;
          if (verbosity > 0)
            solver->report();
          nwork = solver->end_iteration(x, g);
        }
        molpro::cout << "Error={ ";
        for (const auto& e : solver->errors())
          molpro::cout << e << " ";
        molpro::cout << "} after " << solver->statistics().iterations << " iterations, "
                     << solver->statistics().r_creations << " R vectors" << std::endl;
        EXPECT_THAT(solver->errors(), ::testing::Pointwise(::testing::DoubleNear(2 * solver->convergence_threshold()),
                                                           std::vector<double>(nroot, double(0))));
        EXPECT_NEAR(solver->value(), 0, 2e-9);
        const auto nR_creations = solver->statistics().r_creations;
        if (verbosity > 0)
          molpro::cout << "R creations = " << nR_creations << std::endl;
        EXPECT_LE(nR_creations, (nroot + 1) * n_iter);
        std::vector<double> parameters(n), residual(n);
        std::vector<int> roots(1, 0);
        solver->solution(parameters, residual);
        EXPECT_LE(std::sqrt(dot(residual, residual)), solver->convergence_threshold());
        EXPECT_THAT(parameters, ::testing::Pointwise(::testing::DoubleNear(solver->convergence_threshold()),
                                                     std::vector<double>(n, double(1))));
      }
    }
#ifndef NOFORTRAN
    EXPECT_TRUE(test_optimizef(hmat.data(), n) != 0);
#endif
  }
};

TEST_F(OptimizeF, small_quadratic_form) {
  for (int n = 2; n < 31; n+=9) {
    double param = 10;
    load_matrix(n, "", param);
    test_quadratic_form("BFGS", std::to_string(n) + "/" + std::to_string(param));
    test_quadratic_form("SD", std::to_string(n) + "/" + std::to_string(param));
  }
}
auto rosenbrock(const std::vector<double>& x, double a = 1, double b = 100) {
  double f{0};
  std::vector<double> f1(x.size(), 0);
  for (size_t i = 0; i < x.size() - 1; i++) {
    f += std::pow(a - x[i], 2) + b * std::pow(std::pow(x[i], 2) - x[i + 1], 2);
    f1[i] -= 2 * (a - x[i]);
    f1[i] += 4 * b * x[i] * (std::pow(x[i], 2) - x[i + 1]);
    f1[i + 1] -= 2 * b * (std::pow(x[i], 2) - x[i + 1]);
  }
  return std::tuple{x, f, f1};
}
TEST(Optimize, Rosenbrock) {

  for (size_t n = 2; n < 7; n++) {
    auto solver = molpro::linalg::itsolv::create_Optimize<Rvector, Qvector>(
        "BFGS", "convergence_threshold=1e-8,max_size_qspace=6");
    std::vector<double> x(n, -4), g(n);
    x[0] = -3.0;
    for (int iter = 0; iter < 10000; iter++) {
      auto [xx, value, g] = rosenbrock(x);
      //      molpro::cout << "iter=" << iter << ", x=" << x << ", value=" << value << std::endl;
      auto precon = solver->add_vector(x, g, value) > 0;
      //      molpro::cout << "add_vector returns precon=" << precon << ", x=" << x << ", g=" << g[0] << std::endl;
      if (precon)
        for (auto& gg : g)
          gg /= 2;
      if (solver->end_iteration(x, g) == 0)
        break;
      //      molpro::cout << "end_iteration returns x=" << x << ", g=" << g << std::endl;
    }
    molpro::cout << solver->statistics() << std::endl;
    EXPECT_THAT(x, ::testing::Pointwise(::testing::DoubleNear(solver->convergence_threshold()),
                                        std::vector<double>(x.size(), double(1))));
  }
}

TEST(Optimize, trig1d) {

  auto solver =
      molpro::linalg::itsolv::create_Optimize<Rvector, Qvector>("BFGS", "convergence_threshold=1e-8,max_size_qspace=2");
  std::vector<double> x(1), g(1);
//  int nwork = 1;
  x[0] = 1.0;
  for (int iter = 0; iter < 100; iter++) {
    double value = std::sin(x[0]);
    molpro::cout << "iter=" << iter << ", x=" << x << ", value=" << value << std::endl;
    g[0] = std::cos(x[0]);
    auto precon = solver->add_vector(x, g, value) > 0;
    molpro::cout << "add_value returns precon=" << precon << ", x=" << x[0] << ", g=" << g[0] << std::endl;
    //    if (precon)
    //      g[0] = g[0];
    if (solver->end_iteration(x, g) == 0)
      break;
    molpro::cout << "end_iteration returns x=" << x[0] << ", g=" << g[0] << std::endl;
  }
}
