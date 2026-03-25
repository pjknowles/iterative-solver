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
    solver->set_verbosity(molpro::linalg::itsolv::Verbosity::Detailed);
    solver->set_max_iter(100);
    std::vector<Rvector> c = {{1,-1},{-1,1}};
    std::vector<Rvector> g = {{0,0},{0,0}};
    if (not solver->solve(c, g, problem, false))
      std::cout << "failed" << std::endl;
    else
      std::cout << "converged in " << solver->statistics().iterations << " iterations" << std::endl;
    std::vector<int> roots = {0, 1};
    solver->solution(roots, c, g);
  }
  molpro::mpi::finalize();
}

#include <molpro/linalg/itsolv/SolverFactory-implementation.h>
#include <vector>
template class molpro::linalg::itsolv::SolverFactory<std::vector<double>, std::vector<double>,
                                                     std::map<size_t, double>>;
