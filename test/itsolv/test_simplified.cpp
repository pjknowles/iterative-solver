//#include "test.h"
#include "vector_types.h"
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <molpro/iostream.h>
#include <molpro/linalg/itsolv/CastOptions.h>
#include <molpro/linalg/itsolv/SolverFactory.h>
#include <molpro/linalg/itsolv/helper.h>
#include <vector>

// Find lowest eigensolution of a matrix obtained from an external file
// Storage of vectors in-memory via class Rvector

using molpro::linalg::array::Span;
using molpro::linalg::itsolv::CastOptions;
using molpro::linalg::itsolv::CVecRef;
using molpro::linalg::itsolv::cwrap;
using molpro::linalg::itsolv::LinearEigensystem;
using molpro::linalg::itsolv::VecRef;
using molpro::linalg::itsolv::wrap;

struct simplified : ::testing::Test {

  class MyProblem : public molpro::linalg::itsolv::Problem<Rvector> {
  protected:
    //    std::unique_ptr<Qvector> m_diagonals;
    double matrix(int i, int j) const { return i == j ? i + 1 : 0.001 * (i + j); }

  public:
    const size_t n;
    //        MyProblem(const Rvector &diagonals) : m_diagonals(new Qvector(diagonals)) {}
    MyProblem(int n = 10) : n(n) {
      // m_diagonals.reset(new Qvector(n));
    }

    void precondition(const VecRef<Rvector> &action,
                      const std::vector<typename Rvector::value_type> &shift) const override {
      for (size_t k = 0; k < action.size(); k++) {
        auto &a = action[k].get();
        for (size_t i = 0; i < a.size(); i++)
          a[i] /= (matrix(i, i) + shift[k]);
      }
    }

    double residual(const Rvector &v, Rvector &a) const override {
      double value = 0;
      for (size_t i = 0; i < a.size(); i++) {
        a[i] = 0;
        for (size_t j = 0; j < a.size(); j++)
          a[i] += matrix(i, j) * (v[j] - 1);
        value += 0.5 * a[i] * (v[i] - 1);
      }
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
  };
};

TEST_F(simplified, Problem) {
  auto problem = simplified::MyProblem(10);
  auto nonlinear_solver = molpro::linalg::itsolv::create_Optimize<Rvector, Qvector>("BFGS");
  auto linear_solver = molpro::linalg::itsolv::create_LinearEquations<Rvector, Qvector>();
  Rvector v0(problem.n), v1(problem.n);
  EXPECT_TRUE(nonlinear_solver->test_problem(problem,v0,v1,0,1e-9));
  EXPECT_TRUE(linear_solver->test_problem(problem,v0,v1,0,1e-9));
}

TEST_F(simplified, BFGS) {
  auto solver = molpro::linalg::itsolv::create_Optimize<Rvector, Qvector>("BFGS");
  solver->set_convergence_threshold(1e-10);
  auto problem = simplified::MyProblem(10);
  Rvector c(problem.n), g(problem.n);
  EXPECT_TRUE(solver->solve(c, g, problem));
  solver->solution(c, g);
  EXPECT_NEAR(solver->value(), 0, 1e-20);
  EXPECT_THAT(c, ::testing::Pointwise(::testing::DoubleNear(1e-10), Rvector(problem.n, 1)));
  EXPECT_THAT(g, ::testing::Pointwise(::testing::DoubleNear(1e-10), Rvector(problem.n, 0)));
}

TEST_F(simplified, DIIS) {
  auto solver = molpro::linalg::itsolv::create_NonLinearEquations<Rvector, Qvector>("DIIS", "max_size_qspace=4");

  solver->set_convergence_threshold(1e-12);
  auto problem = simplified::MyProblem(10);
  Rvector c(problem.n), g(problem.n);
  EXPECT_TRUE(solver->solve(c, g, problem));
  solver->solution(c, g);
  EXPECT_THAT(c, ::testing::Pointwise(::testing::DoubleNear(1e-10), Rvector(problem.n, 1)));
  EXPECT_THAT(g, ::testing::Pointwise(::testing::DoubleNear(solver->convergence_threshold()), Rvector(problem.n, 0)));
}

TEST_F(simplified, LinearEigensystem) {
  auto solver = molpro::linalg::itsolv::create_LinearEigensystem<Rvector, Qvector>("Davidson");
  solver->set_convergence_threshold(1e-10);
  auto problem = simplified::MyProblem(10);
  Rvector c(problem.n, 0), g(problem.n);
  c[0] = 1;
  EXPECT_TRUE(solver->solve(c, g, problem));
  solver->solution(c, g);
  EXPECT_THAT(g, ::testing::Pointwise(::testing::DoubleNear(1e-10), Rvector(problem.n, 0)));
}

template <typename s>
std::ostream &operator<<(std::ostream &o, const std::vector<s> &v) {
  for (const auto &e : v)
    o << " " << e;
  return o;
}
TEST_F(simplified, LinearEquations) {
  auto solver = molpro::linalg::itsolv::create_LinearEquations<Rvector, Qvector>("Davidson");
  solver->set_convergence_threshold(1e-10);
  auto problem = simplified::MyProblem(10);
  Rvector c(problem.n, 1), g(problem.n), rhs(problem.n);
  problem.action(molpro::linalg::itsolv::cwrap_arg(c), molpro::linalg::itsolv::wrap_arg(rhs));
  solver->add_equations(molpro::linalg::itsolv::cwrap_arg(rhs));
  c = rhs;
  EXPECT_TRUE(solver->solve(c, g, problem));
  solver->solution(c, g);
  EXPECT_THAT(c, ::testing::Pointwise(::testing::DoubleNear(1e-10), Rvector(problem.n, 1)));
  EXPECT_THAT(g, ::testing::Pointwise(::testing::DoubleNear(1e-10), Rvector(problem.n, 0)));
  problem.action(molpro::linalg::itsolv::cwrap_arg(c), molpro::linalg::itsolv::wrap_arg(g));
}