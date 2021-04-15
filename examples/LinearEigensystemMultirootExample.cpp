#include "ExampleProblem.h"
#include <iostream>
#include <molpro/linalg/itsolv/SolverFactory.h>

int main(int argc, char* argv[]) {
  auto problem = ExampleProblem(argc > 1 ? std::stoi(argv[1]) : 20);
  using Rvector = ExampleProblem::container_t;
  using Qvector = Rvector;
  auto solver = molpro::linalg::itsolv::create_LinearEigensystem<Rvector, Qvector>("Davidson");
  solver->set_n_roots(argc > 2 ? std::stoi(argv[2]) : 2);
  solver->set_verbosity(molpro::linalg::itsolv::Verbosity::Summary);
  solver->set_max_iter(100);
  std::vector<Rvector> c, g;
  auto nbuffer = argc > 3 ? std::stoi(argv[3]) : solver->n_roots();
  for (int root = 0; root < nbuffer; root++) {
    c.emplace_back(problem.n, 0);
    c.back()[root] = 1;
    g.emplace_back(problem.n);
  }
  if (not solver->solve(c, g, problem))
    std::cout << "failed" << std::endl;
  else
    std::cout << "converged in " << solver->statistics().iterations << " iterations" << std::endl;
  std::vector<int> roots;
  std::iota(roots.begin(), roots.end(), 0);
  solver->solution(roots, c, g);
  for (const auto& ev : solver->eigenvalues())
    std::cout << "Final eigenvalue: " << ev << std::endl;
}