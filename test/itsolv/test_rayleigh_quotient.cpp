#include "vector_types.h"
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <molpro/iostream.h>
#include <molpro/linalg/itsolv/CastOptions.h>
#include <molpro/linalg/itsolv/SolverFactory.h>
#include <molpro/linalg/itsolv/helper.h>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

using molpro::linalg::array::Span;
using molpro::linalg::itsolv::CastOptions;
using molpro::linalg::itsolv::CVecRef;
using molpro::linalg::itsolv::cwrap;
using molpro::linalg::itsolv::LinearEigensystem;
using molpro::linalg::itsolv::VecRef;
using molpro::linalg::itsolv::wrap;

template <typename s>
std::ostream &operator<<(std::ostream &o, const std::vector<s> &v) {
  for (const auto &e : v)
    o << " " << e;
  return o;
}

struct RayleighQuotient : ::testing::Test {

  class MyProblem : public molpro::linalg::itsolv::Problem<Rvector> {
  protected:
    //    double matrix(int i, int j) const { return i == j ? i + 1 + rho : rho; }
    Eigen::MatrixXd matrix;

  public:
    const size_t n;
    const double rho;
    MyProblem(int n = 10, double rho = 0.1) : n(n), rho(rho) {
      matrix.resize(n,n);
      for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
          matrix(i, j) = i == j ? i + 1 + rho : rho;
    }

    bool diagonals(container_t &d) const override {
      for (size_t i = 0; i < d.size(); i++)
        d[i] = matrix(i, i);
      return true;
    }

    double residual(const Rvector &v, Rvector &a) const override {
      double norm = std::inner_product(v.begin(), v.end(), v.begin(), double(0));
      for (size_t i = 0; i < a.size(); i++) {
        a[i] = 0;
        for (size_t j = 0; j < a.size(); j++)
          a[i] += matrix(i, j) * v[j];
      }
      double value = std::inner_product(v.begin(), v.end(), a.begin(), double(0)) / norm;
      for (size_t i = 0; i < a.size(); i++)
        a[i] = 2 * (a[i] - value * v[i]) / norm;
//      std::cout << "residual v "<<v<<std::endl;
//      std::cout << "residual a "<<a<<std::endl;
//      std::cout << "residual value "<<value<<std::endl;
      return value;
    }

    void action(const CVecRef<Rvector> &parameters, const VecRef<Rvector> &actions) const override {
      for (size_t k = 0; k < parameters.size(); k++) {
        const auto &v = parameters[k].get();
        auto &a = actions[k].get();
        for (size_t i = 0; i < a.size(); i++) {
          a[i] = 0;
          for (size_t j = 0; j < a.size(); j++)
            a[i] += matrix(i, j) * v[j];
        }
      }
    }

    bool test_parameters(unsigned int instance, Rvector &parameters) const override {
      //      std::cout << "test_parameters "<<instance<<" "<<parameters.size()<<std::endl;
      parameters.assign(parameters.size(), 1.0);
      if (instance == 0)
        return true;
      if (instance <= parameters.size()) {
        parameters[instance - 1] += 0.0001;
        //        std::cout << "set parameters "<<parameters[0]<<","<<parameters[1]<<",..."<<std::endl;
        return true;
      }
      return false;
    }

    Eigen::MatrixXd eigenvectors() const {
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(matrix);
      return solver.eigenvectors();
    }

    Eigen::MatrixXd eigenvalues() const {
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(matrix);
      return solver.eigenvalues();
    }
    std::vector<double> lowest_eigenvector() const {
      auto eigenvector = eigenvectors().col(0).eval();
      std::vector<double> expected(eigenvector.data(), eigenvector.data()+n);
      return expected;
    }
  };
};

TEST_F(RayleighQuotient, Problem) {
  auto problem = RayleighQuotient::MyProblem(4, 0.01);
  auto nonlinear_solver = molpro::linalg::itsolv::create_Optimize<Rvector, Qvector>("BFGS");
  auto linear_solver = molpro::linalg::itsolv::create_LinearEquations<Rvector, Qvector>();
  Rvector v0(problem.n), v1(problem.n);
  EXPECT_TRUE(nonlinear_solver->test_problem(problem, v0, v1, 0, 1e-9));
  EXPECT_TRUE(linear_solver->test_problem(problem, v0, v1, 0, 1e-9));
}

TEST_F(RayleighQuotient, BFGS) {
  auto solver = molpro::linalg::itsolv::create_Optimize<Rvector, Qvector>("BFGS");
  solver->set_convergence_threshold(1e-10);
  auto problem = RayleighQuotient::MyProblem(4, 0.01);
  Rvector c(problem.n), g(problem.n);
  c.assign(problem.n, double(10));
  EXPECT_TRUE(solver->solve(c, g, problem));
  solver->solution(c, g);
  auto norm = std::sqrt(std::inner_product(c.begin(), c.end(), c.begin(), double(0)));
  std::transform(c.begin(), c.end(), c.begin(), [norm](const double a) { return a / norm; });
  EXPECT_NEAR(solver->value(), problem.eigenvalues()(0), 1e-14);
  EXPECT_THAT(c, ::testing::Pointwise(::testing::DoubleNear(1e-8), problem.lowest_eigenvector()));
  EXPECT_THAT(g, ::testing::Pointwise(::testing::DoubleNear(solver->convergence_threshold()*10), Rvector(problem.n, 0)));
}

TEST_F(RayleighQuotient, DIIS) {
  auto solver = molpro::linalg::itsolv::create_NonLinearEquations<Rvector, Qvector>("DIIS", "max_size_qspace=4");

  solver->set_convergence_threshold(1e-12);
  auto problem = RayleighQuotient::MyProblem(4, 0.01);
  Rvector c(problem.n), g(problem.n);
  c.assign(problem.n, double(0));
  c[0]=1;
  EXPECT_TRUE(solver->solve(c, g, problem));
  solver->solution(c, g);
  auto norm = std::sqrt(std::inner_product(c.begin(), c.end(), c.begin(), double(0)));
  std::transform(c.begin(), c.end(), c.begin(), [norm](const double a) { return a / norm; });
  EXPECT_THAT(c, ::testing::Pointwise(::testing::DoubleNear(1e-9), problem.lowest_eigenvector()));
  EXPECT_THAT(g, ::testing::Pointwise(::testing::DoubleNear(solver->convergence_threshold()*10), Rvector(problem.n, 0)));
}

TEST_F(RayleighQuotient, LinearEigensystem) {
  auto solver = molpro::linalg::itsolv::create_LinearEigensystem<Rvector, Qvector>("Davidson");
  solver->set_convergence_threshold(1e-10);
  auto problem = RayleighQuotient::MyProblem(4, 0.01);
  Rvector c(problem.n, 0), g(problem.n);
  c[0] = 1;
  EXPECT_TRUE(solver->solve(c, g, problem));
  solver->solution(c, g);
  EXPECT_NEAR(solver->eigenvalues()[0], problem.eigenvalues()(0), 1e-10);
  auto norm = std::sqrt(std::inner_product(c.begin(), c.end(), c.begin(), double(0)));
  std::transform(c.begin(), c.end(), c.begin(), [norm](const double a) { return a / norm; });
  EXPECT_THAT(c, ::testing::Pointwise(::testing::DoubleNear(1e-9), problem.lowest_eigenvector()));
  EXPECT_THAT(g, ::testing::Pointwise(::testing::DoubleNear(solver->convergence_threshold()*10), Rvector(problem.n, 0)));
}


TEST_F(RayleighQuotient, LinearEquations) {
  auto solver = molpro::linalg::itsolv::create_LinearEquations<Rvector, Qvector>("Davidson");
  solver->set_convergence_threshold(1e-10);
  auto problem = RayleighQuotient::MyProblem(4, 0.01);
  Rvector c(problem.n, 1), g(problem.n);
  auto rhs = problem.lowest_eigenvector();
  auto eval = problem.eigenvalues()(0);
  std::transform(rhs.begin(),rhs.end(),rhs.begin(),[eval](double a){return a*eval;});
  solver->add_equations(molpro::linalg::itsolv::cwrap_arg(rhs));
//  c = rhs;
  EXPECT_TRUE(solver->solve(c, g, problem));
  solver->solution(c, g);
  auto norm = std::sqrt(std::inner_product(c.begin(), c.end(), c.begin(), double(0)));
  std::transform(c.begin(), c.end(), c.begin(), [norm](const double a) { return a / norm; });
  EXPECT_THAT(c, ::testing::Pointwise(::testing::DoubleNear(1e-9), problem.lowest_eigenvector()));
  EXPECT_THAT(g, ::testing::Pointwise(::testing::DoubleNear(solver->convergence_threshold()*10), Rvector(problem.n, 0)));
}