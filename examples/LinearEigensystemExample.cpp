#include "ExampleProblem.h"
#include <iostream>
#include <molpro/linalg/itsolv/SolverFactory.h>
#include <molpro/mpi.h>

int main(int argc, char* argv[]) {
  molpro::mpi::init();
  {
    auto problem = ExampleProblem(argc > 1 ? std::stoi(argv[1]) : 20);
    using Rvector = ExampleProblem::container_t;
    auto solver = molpro::linalg::itsolv::create_LinearEigensystem<Rvector>("Davidson");
    solver->set_n_roots(argc > 2 ? std::stoi(argv[2]) : 2);
    //  solver->set_verbosity(molpro::linalg::itsolv::Verbosity::Detailed);
    solver->set_max_iter(100);
    Rvector c(problem.n), g(problem.n);
    if (not solver->solve(c, g, problem, true))
      std::cout << "failed" << std::endl;
    else
      std::cout << "converged in " << solver->statistics().iterations << " iterations" << std::endl;
    solver->solution(c, g);
    for (const auto& ev : solver->eigenvalues())
      std::cout << "Final eigenvalue: " << ev << std::endl;
  }
  molpro::mpi::finalize();
}

#include <molpro/linalg/itsolv/SolverFactory-implementation.h>
#include <vector>
template class molpro::linalg::itsolv::SolverFactory<std::vector<double>, std::vector<double>,
                                                     std::map<size_t, double>>;
