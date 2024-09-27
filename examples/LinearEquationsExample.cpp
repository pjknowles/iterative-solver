#include "ExampleProblem.h"
#include <iostream>
#include <molpro/linalg/itsolv/SolverFactory.h>
#include <molpro/mpi.h>

int main(int argc, char* argv[]) {
  molpro::mpi::init();
  {
    auto problem = ExampleProblem(argc > 1 ? std::stoi(argv[1]) : 20);
    using Rvector = ExampleProblem::container_t;

    // old
    auto solver = molpro::linalg::itsolv::create_LinearEquations<Rvector>("Davidson");
    solver->set_verbosity(molpro::linalg::itsolv::Verbosity::Summary);
    Rvector c(problem.n), g(problem.n);
    if (not solver->solve(c, g, problem, true))
      std::cout << "failed" << std::endl;
    else
      std::cout << "converged in " << solver->statistics().iterations << " iterations" << std::endl;
    solver->solution(c, g);

    // new
    auto [success, solver2] = molpro::linalg::itsolv::Solve_LinearEquations<Rvector>(c, g, problem, 1);
    // including default arguments:
    //    auto [success, solver2] = molpro::linalg::itsolv::Solve_LinearEquations<Rvector>(c, g, problem, 1, true,
    //    "Davidson", "");
    if (!success)
      std::cout << "failed" << std::endl;
    else
      std::cout << "converged in " << solver->statistics().iterations << " iterations" << std::endl;
    std::cout << solver2->statistics() << std::endl;
    std::cout << "solution";
    for (const auto& e : c)
      std::cout << " " << e;
    std::cout << std::endl;
    std::cout << "residual";
    for (const auto& e : g)
      std::cout << " " << e;
    std::cout << std::endl;
  }
  molpro::mpi::finalize();
}

#include <molpro/linalg/itsolv/SolverFactory-implementation.h>
#include <vector>
template class molpro::linalg::itsolv::SolverFactory<std::vector<double>, std::vector<double>,
                                                     std::map<size_t, double>>;
