#include "ExampleProblem.h"
#include <iostream>
#include <molpro/linalg/itsolv/SolverFactory.h>
#include <molpro/mpi.h>

int main(int argc, char* argv[]) {
  molpro::mpi::init();
  {
    auto problem = ExampleProblem(argc > 1 ? std::stoi(argv[1]) : 20);
    using Rvector = ExampleProblem::container_t;
    auto solver = molpro::linalg::itsolv::create_LinearEquations<Rvector>("Davidson");
    solver->set_verbosity(molpro::linalg::itsolv::Verbosity::Summary);
    Rvector c(problem.n), g(problem.n), rhs(problem.n, 1);
    solver->add_equations(molpro::linalg::itsolv::cwrap_arg(rhs));
    if (not solver->solve(c, g, problem, true))
      std::cout << "failed" << std::endl;
    else
      std::cout << "converged in " << solver->statistics().iterations << " iterations" << std::endl;
    solver->solution(c, g);
  }
  molpro::mpi::finalize();
}

#include <molpro/linalg/itsolv/SolverFactory-implementation.h>
#include <vector>
template class molpro::linalg::itsolv::SolverFactory<std::vector<double>, std::vector<double>,
                                                     std::map<size_t, double>>;
