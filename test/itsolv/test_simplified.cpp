//#include "test.h"
#include "vector_types.h"
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <molpro/iostream.h>
#include <molpro/linalg/itsolv/CastOptions.h>
#include <molpro/linalg/itsolv/SolverFactory.h>
#include <molpro/linalg/itsolv/helper.h>
#include <algorithm>
#include <numeric>
#include <vector>

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

    bool diagonals(container_t &d) const override {
      for (size_t i = 0; i < d.size(); i++)
        d[i] = matrix(i,i);
    return true;
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

  class PSpaceProblem : public MyProblem {
  public:
    using MyProblem::MyProblem;

    std::vector<double> pp_action_matrix(const std::vector<Pvector> &pparams) const override {
      std::vector<double> result;
      for (const auto &pi : pparams)
        for (const auto &pj : pparams) {
          result.push_back(0);
          for (const auto &pie : pi)
            for (const auto &pje : pj)
              result.back() += matrix(pie.first, pje.first) * pie.second * pje.second;
        }
      return result;
    }

    void p_action(const std::vector<std::vector<double>> &p_coefficients, const CVecRef<Pvector> &pparams,
                  const VecRef<Rvector> &actions) const override {
      for (size_t k = 0; k < p_coefficients.size(); k++) {
        auto &a = actions[k].get();
        for (size_t pi = 0; pi < pparams.size(); pi++)
          for (const auto &pie : pparams[pi].get())
            for (size_t i = 0; i < a.size(); i++)
              a[i] += matrix(i, pie.first) * pie.second * p_coefficients[k][pi];
      }
    }
  };

  class CountingProblem : public MyProblem {
  public:
    using MyProblem::MyProblem;
    mutable std::vector<size_t> action_sizes; //!< number of vectors passed to each call of action()

    void action(const CVecRef<Rvector> &parameters, const VecRef<Rvector> &actions) const override {
      action_sizes.push_back(parameters.size());
      MyProblem::action(parameters, actions);
    }
  };

  class RecordingLogger : public molpro::linalg::itsolv::Logger {
  public:
    mutable std::vector<std::string> messages;

    size_t count(std::string_view text) const {
      return std::count_if(messages.begin(), messages.end(),
                           [text](const std::string &m) { return m.find(text) != std::string::npos; });
    }

  protected:
    void default_message_handler(std::string_view, molpro::linalg::itsolv::log::Severity,
                                 molpro::linalg::itsolv::log::Verbosity, std::string_view msg) const override {
      messages.emplace_back(msg);
    }
  };

  // solves for the lowest nroots eigenpairs and checks that they converged
  template <class Problem>
  static void solve_eigensystem(const Problem &problem, size_t nroots, size_t nbuffers, int max_p,
                                double p_threshold = 1e10) {
    auto solver = molpro::linalg::itsolv::create_LinearEigensystem<Rvector, Qvector>("Davidson");
    solver->set_hermiticity(true);
    solver->set_n_roots(nroots);
    solver->set_max_p(max_p);
    solver->set_p_threshold(p_threshold);
    std::vector<Rvector> c(nbuffers, Rvector(problem.n)), g(nbuffers, Rvector(problem.n));
    ASSERT_TRUE(solver->solve(c, g, problem, true));
    ASSERT_EQ(solver->eigenvalues().size(), nroots);
    std::vector<int> roots(nroots);
    std::iota(roots.begin(), roots.end(), 0);
    c.resize(nroots, Rvector(problem.n));
    g.resize(nroots, Rvector(problem.n));
    solver->solution(roots, c, g);
    for (size_t i = 0; i < nroots; i++)
      EXPECT_THAT(g[i], ::testing::Pointwise(::testing::DoubleNear(solver->convergence_threshold()),
                                             Rvector(problem.n, 0)));
  }
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
TEST_F(simplified, LinearEigensystem_max_iter_zero) {
  auto solver = molpro::linalg::itsolv::create_LinearEigensystem<Rvector, Qvector>("Davidson");
  solver->set_max_iter(0);
  auto problem = simplified::MyProblem(10);
  Rvector c(problem.n, 0), g(problem.n);
  c[0] = 1;
  EXPECT_FALSE(solver->solve(c, g, problem));
  EXPECT_TRUE(solver->errors().empty());
}

TEST_F(simplified, LinearEigensystem_pspace) {
  auto problem = simplified::PSpaceProblem(10);
  solve_eigensystem(problem, 3, 3, 4);
}

TEST_F(simplified, LinearEigensystem_pspace_max_p_smaller_than_n_roots) {
  auto problem = simplified::PSpaceProblem(10);
  solve_eigensystem(problem, 3, 3, 2);
}

TEST_F(simplified, LinearEigensystem_pspace_threshold_smaller_than_n_roots) {
  auto problem = simplified::PSpaceProblem(10);
  solve_eigensystem(problem, 3, 3, 5, 0.5);
}

TEST_F(simplified, LinearEigensystem_pspace_too_small_zero_parameters) {
  // the parameters are not used with max_p > 0, so they need not hold an initial guess
  const size_t nroots = 3;
  auto problem = simplified::PSpaceProblem(10);
  auto solver = molpro::linalg::itsolv::create_LinearEigensystem<Rvector, Qvector>("Davidson");
  solver->set_hermiticity(true);
  solver->set_n_roots(nroots);
  solver->set_max_p(2);
  std::vector<Rvector> c(nroots, Rvector(problem.n, 0)), g(nroots, Rvector(problem.n));
  EXPECT_TRUE(solver->solve(c, g, problem));
  const auto eigenvalues = solver->eigenvalues();
  ASSERT_EQ(eigenvalues.size(), nroots);
  for (size_t i = 0; i < nroots; i++)
    EXPECT_NEAR(eigenvalues[i], i + 1, 1e-3);
}

TEST_F(simplified, LinearEigensystem_surplus_buffers) {
  const size_t nroots = 3;
  auto problem = simplified::CountingProblem(10);
  auto solver = molpro::linalg::itsolv::create_LinearEigensystem<Rvector, Qvector>("Davidson");
  solver->set_n_roots(nroots);
  std::vector<Rvector> c(2 * nroots, Rvector(problem.n)), g(2 * nroots, Rvector(problem.n));
  EXPECT_TRUE(solver->solve(c, g, problem, true));
  ASSERT_FALSE(problem.action_sizes.empty());
  EXPECT_EQ(problem.action_sizes.front(), nroots);
  // every computed action enters the subspace
  EXPECT_EQ(std::accumulate(problem.action_sizes.begin(), problem.action_sizes.end(), size_t{0}),
            solver->statistics().r_creations);
  const auto eigenvalues = solver->eigenvalues();
  ASSERT_EQ(eigenvalues.size(), nroots);
  for (size_t i = 0; i < nroots; i++)
    EXPECT_NEAR(eigenvalues[i], i + 1, 1e-3);
}

TEST_F(simplified, LinearEigensystem_solve_converged) {
  auto problem = simplified::CountingProblem(10);
  auto solver = molpro::linalg::itsolv::create_LinearEigensystem<Rvector, Qvector>("Davidson");
  solver->set_n_roots(2);
  std::vector<Rvector> c(2, Rvector(problem.n)), g(2, Rvector(problem.n));
  ASSERT_TRUE(solver->solve(c, g, problem, true));
  ASSERT_TRUE(solver->working_set().empty());
  problem.action_sizes.clear();
  EXPECT_TRUE(solver->solve(c, g, problem));
  // nothing is left to compute, so action() is not called
  EXPECT_TRUE(problem.action_sizes.empty());
  EXPECT_TRUE(solver->working_set().empty());
}

TEST_F(simplified, LinearEigensystem_pspace_no_initial_guess) {
  // the P space replaces the initial guess, so none is generated
  auto problem = simplified::PSpaceProblem(10);
  auto solver = molpro::linalg::itsolv::create_LinearEigensystem<Rvector, Qvector>("Davidson");
  auto logger = std::make_shared<RecordingLogger>();
  solver->set_logger(logger);
  solver->set_hermiticity(true);
  solver->set_n_roots(2);
  solver->set_max_p(4);
  std::vector<Rvector> c(2, Rvector(problem.n)), g(2, Rvector(problem.n));
  EXPECT_TRUE(solver->solve(c, g, problem, true));
  EXPECT_EQ(logger->count("P-space dimension"), 1);
  EXPECT_EQ(logger->count("Initial guess generated"), 0);
}
