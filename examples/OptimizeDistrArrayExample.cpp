#include "ExampleProblemDistrArray.h"
#include <iostream>
#include <molpro/linalg/array/DistrArrayFile.h>
#include <molpro/linalg/array/DistrArraySpan.h>
#include <molpro/linalg/array/util/Distribution.h>
#include <molpro/linalg/itsolv/SolverFactory.h>
#include <molpro/mpi.h>

int main(int argc, char* argv[]) {
  molpro::mpi::init();
  using Rvector = molpro::linalg::array::DistrArraySpan;
  auto solver = molpro::linalg::itsolv::create_Optimize<Rvector, molpro::linalg::array::DistrArrayFile>("BFGS");
  auto problem = ExampleProblem(argc > 1 ? std::stoi(argv[1]) : 20);
  solver->set_verbosity(molpro::mpi::rank_global() == 0 ? molpro::linalg::itsolv::Verbosity::Summary
                                                        : molpro::linalg::itsolv::Verbosity::None);
  solver->set_max_iter(100);
  auto distribution = molpro::linalg::array::util::make_distribution_spread_remainder<Rvector::index_type>(
      problem.n, molpro::mpi::size_global());
  auto local_range = distribution.range(molpro::mpi::rank_global());
  problem.m_vectors.front().assign(problem.n, 0);
  problem.m_vectors.front()[0] = 1;
  std::vector<Rvector::value_type> action_vector(local_range.second - local_range.first);
  Rvector c(std::make_unique<molpro::linalg::array::util::Distribution<Rvector::index_type>>(distribution),
                 molpro::linalg::array::Span<Rvector::value_type>(problem.m_vectors.front().data() + local_range.first,
                                                                  local_range.second - local_range.first));
  Rvector g(std::make_unique<molpro::linalg::array::util::Distribution<Rvector::index_type>>(distribution),
                 molpro::linalg::array::Span<Rvector::value_type>(action_vector.data(), action_vector.size()));
  if (not solver->solve(c, g, problem))
    std::cout << "failed" << std::endl;
  else if (molpro::mpi::rank_global() == 0) {
    std::cout << "converged in " << solver->statistics().iterations << " iterations" << std::endl;
  }
  molpro::mpi::finalize();
}

#include <molpro/linalg/itsolv/SolverFactory-implementation.h>
template class molpro::linalg::itsolv::SolverFactory<molpro::linalg::array::DistrArraySpan,
                                                     molpro::linalg::array::DistrArrayFile, std::map<size_t, double>>;
