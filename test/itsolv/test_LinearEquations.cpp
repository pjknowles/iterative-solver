#include "test.h"

#include <Eigen/Eigenvalues>
#include <numeric>

#include "vector_types.h"
#include <molpro/linalg/itsolv/CastOptions.h>
#include <molpro/linalg/itsolv/SolverFactory.h>
#include <molpro/linalg/itsolv/helper.h>

using molpro::linalg::itsolv::CastOptions;
using molpro::linalg::itsolv::LinearEquations;

class Problem_ : public molpro::linalg::itsolv::Problem<Rvector> {
public:
  Problem_(int n, int nroot = 1) {
    this->n = n;
    matrix.resize(n, n);
    rhs.resize(n, nroot);
    expected_solution.resize(n, nroot);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j)
        matrix(i, j) = i + j + 1;
      matrix(i, i) += 1;
    }
    for (int i = 0; i < nroot; ++i)
      for (int j = 0; j < n; ++j) {
        rhs(j, i) = (i + 1) * (n * (n + 1) / 2 + j * n + 1);
        expected_solution(j, i) = i + 1;
      }
  }
  int n;
  Eigen::MatrixXd matrix;
  Eigen::MatrixXd rhs;
  Eigen::MatrixXd expected_solution;
  void action(const CVecRef<Rvector>& parameters, const VecRef<Rvector>& action) const override {
    for (int v = 0; v <parameters.size(); ++v) {
      Eigen::Map<Eigen::VectorXd>(action[v].get().data(), n) =
          matrix * Eigen::Map<const Eigen::VectorXd>(parameters[v].get().data(), n);
//      std::cout << "parameters: " << Eigen::Map<const Eigen::VectorXd>(parameters[v].get().data(), n) << std::endl;
//      std::cout << "action: " << Eigen::Map<const Eigen::VectorXd>(action[v].get().data(), n) << std::endl;
    }
  }
  bool diagonals(Rvector& d) const override {
    for (int i = 0; i < n; ++i)
      d[i] = matrix(i, i);
    return true;
  }

  bool RHS(Rvector& r, unsigned int instance) const override {
    if (instance < rhs.cols()) {
      std::copy(rhs.col(instance).begin(), rhs.col(instance).end(), r.begin());
      return true;
    } else
      return false;
  }
};

#ifndef NOFORTRAN
extern "C" int test_linearequationsf(const double* matrix, const double* rhs, size_t n, size_t np, size_t nroot,
                                     int hermitian, double augmented_hessian);
#endif

TEST(LinearEquation, symmetric_system) {

  size_t n_min = 3, n_max = 33;
  double augmented_hessian = 0;
  for (int n = n_min; n <= n_max; n += 3) {
    for (int nroot = 1; nroot <= n && nroot < 14; nroot++) {
      std::cout << "n=" << n << ", nroot=" << nroot << std::endl;
      Problem_ problem(n, nroot);
//      std::cout << "rhs\n" << problem.rhs << std::endl;
//      std::cout << "expected_solution\n" << problem.expected_solution << std::endl;
//      std::cout << "recreated rhs\n" <<  problem.matrix*problem.expected_solution  << std::endl;
      ASSERT_EQ(  problem.matrix*problem.expected_solution, problem.rhs);
      auto solver = molpro::linalg::itsolv::create_LinearEquations<Rvector, Qvector, Pvector>();
      for (int root = 0; false and root < nroot; ++root) {
        std::vector<double> r(problem.n);
        for (int i = 0; i < problem.n; ++i)
          r[i] = problem.rhs(root, i);
        solver->add_equations(r);
      }
      std::vector<std::vector<double>> parameters(nroot);
      for (auto& v : parameters)
        v.resize(n);
      std::vector<std::vector<double>> actions(nroot);
      for (auto& v : actions)
        v.resize(n);
      solver->set_convergence_threshold(1e-10);
      solver->solve(parameters, actions, problem, true);
      std::vector<int> roots(nroot);
      std::iota(roots.begin(), roots.end(), 0);
      solver->solution(roots, parameters, actions);
      for (int root = 0; root < nroot; ++root) {
//        std::cout << "Solution\n" << Eigen::Map<const Eigen::VectorXd>(parameters[root].data(), n) << std::endl;
        for (int i = 0; i < n; ++i)
          EXPECT_NEAR(parameters[root][i], problem.expected_solution(i, root), 1e-5);
      }
      //            auto options = set_options(solver, mat.rows(), nroot, 0, hermitian, augmented_hessian);
#ifndef NOFORTRAN
      EXPECT_TRUE(test_linearequationsf(problem.matrix.data(), problem.rhs.data(), n, 0, nroot, true, augmented_hessian) != 0);
#endif
    }
  }
}