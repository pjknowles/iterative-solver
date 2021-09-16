#include "ExampleProblemDistrArray.h"
#include <iostream>
#include <molpro/linalg/array/DistrArrayFile.h>
#include <molpro/linalg/array/DistrArraySpan.h>
#include <molpro/linalg/array/util/Distribution.h>
#include <molpro/linalg/itsolv/SolverFactory.h>
#include <molpro/mpi.h>

int main(int argc, char* argv[]) {
  molpro::mpi::init();
  {
    using Rvector = molpro::linalg::array::DistrArraySpan;
    auto solver =
        molpro::linalg::itsolv::create_LinearEigensystem<Rvector, molpro::linalg::array::DistrArrayFile>("Davidson");
    solver->set_n_roots(argc > 2 ? std::stoi(argv[2]) : 2);
    auto problem = ExampleProblem(argc > 1 ? std::stoi(argv[1]) : 50, argc > 3 ? std::stoi(argv[3]) : solver->n_roots());
    solver->set_verbosity(molpro::mpi::rank_global() == 0 ? molpro::linalg::itsolv::Verbosity::Summary
                                                          : molpro::linalg::itsolv::Verbosity::None);
    //  solver->set_verbosity(molpro::linalg::itsolv::Verbosity::Detailed);
    solver->set_max_iter(100);
    solver->set_max_p(5);
    //  solver->set_p_threshold(3.7);
    solver->set_hermiticity(true);
    auto distribution = molpro::linalg::array::util::make_distribution_spread_remainder<Rvector::index_type>(
        problem.n, molpro::mpi::size_global());
    auto local_range = distribution.range(molpro::mpi::rank_global());
    std::vector<std::vector<Rvector::value_type>> action_vectors;
    std::vector<Rvector> c, g;
    for (size_t i = 0; i < problem.m_vectors.size(); i++) {
      action_vectors.emplace_back(local_range.second - local_range.first);
      c.emplace_back(std::make_unique<molpro::linalg::array::util::Distribution<Rvector::index_type>>(distribution),
                     molpro::linalg::array::Span<Rvector::value_type>(problem.m_vectors[i].data() + local_range.first,
                                                                      local_range.second - local_range.first));
      g.emplace_back(std::make_unique<molpro::linalg::array::util::Distribution<Rvector::index_type>>(distribution),
                     molpro::linalg::array::Span<Rvector::value_type>(action_vectors.back().data(),
                                                                      local_range.second - local_range.first));
    }
    if (not solver->solve(c, g, problem, true))
      std::cout << "failed" << std::endl;
    else if (molpro::mpi::rank_global() == 0) {
      std::cout << "converged in " << solver->statistics().iterations << " iterations" << std::endl;
    }
  }
  molpro::mpi::finalize();
}

#include <molpro/linalg/itsolv/SolverFactory-implementation.h>
template class molpro::linalg::itsolv::SolverFactory<molpro::linalg::array::DistrArraySpan,
                                                     molpro::linalg::array::DistrArrayFile, std::map<size_t, double>>;
