#include "ExampleProblem.h"
//#include <molpro/linalg/array/DistrArrayFile.h>
#include <iostream>
#include <molpro/linalg/itsolv/SolverFactory.h>

int main(int argc, char* argv[]) {
  auto problem = ExampleProblem(argc > 1 ? std::stoi(argv[1]) : 20);
  using Rvector = ExampleProblem::container_t;
  using Qvector = ExampleProblem::container_t;
  // using Qvector = molpro::linalg::array::DistrArrayFile;
  auto solver = molpro::linalg::itsolv::create_LinearEigensystem<Rvector, Qvector>("Davidson");
  solver->set_n_roots(argc > 2 ? std::stoi(argv[2]) : 2);
  solver->set_verbosity(molpro::linalg::itsolv::Verbosity::Summary);
  solver->set_max_iter(100);
  Rvector c(problem.n, 0), g(problem.n);
  c[0] = 1;
  if (not solver->solve(c, g, problem))
    std::cout << "failed" << std::endl;
  else
    std::cout << "converged in " << solver->statistics().iterations << " iterations" << std::endl;
  solver->solution(c, g);
  for (const auto& ev : solver->eigenvalues())
    std::cout << "Final eigenvalue: " << ev << std::endl;
}