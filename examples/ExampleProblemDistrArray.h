#ifndef LINEARALGEBRA_EXAMPLES_EXAMPLEPROBLEM_H_
#define LINEARALGEBRA_EXAMPLES_EXAMPLEPROBLEM_H_
#include <molpro/linalg/array/DistrArraySpan.h>
#include <molpro/linalg/array/util/Distribution.h>
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/mpi.h>
#include <vector>

class ExampleProblem : public molpro::linalg::itsolv::Problem<molpro::linalg::array::DistrArraySpan> {

protected:
  double matrix(int i, int j) const { return i == j ? i + 1 : 0.001 * ((i + j) % n); }

public:
  using Problem::container_t;
  using Problem::value_t;
  const size_t n;
  std::vector<std::vector<double>> m_vectors;
  ExampleProblem(int n = 10, size_t nvector = 1) : n(n) {
    for (size_t i = 0; i < nvector; i++) {
      m_vectors.emplace_back(n);
    }
  }

  void precondition(const VecRef<container_t>& action, const std::vector<value_t>& shift) const override {
    for (int k = 0; k < action.size(); k++) {
      auto action_local_buffer = action[k].get().local_buffer();
      auto& distribution = action[k].get().distribution();
      auto range = distribution.range(molpro::mpi::rank_global());
      for (auto i = range.first; i < range.second; i++)
        (*action_local_buffer)[i - range.first] /= (matrix(i, i) - shift[k] + 1e-15);
    }
  }

  value_t residual(const container_t& v, container_t& a) const override {
    auto action_local_buffer_pointer = a.local_buffer();
    auto& action_local_buffer = *action_local_buffer_pointer;
    auto& distribution = a.distribution();
    assert(distribution.compatible(v.distribution()));
    assert(v.local_buffer().get()->data() <= m_vectors.front().data());
    assert(v.local_buffer().get()->data() + n > m_vectors[0].data());
#ifdef HAVE_MPI_H
    for (int rank = 0; rank < molpro::mpi::size_global(); rank++)
      MPI_Bcast((void*)(m_vectors.front().data() + distribution.range(rank).first),
                distribution.range(rank).second - distribution.range(rank).first, MPI_DOUBLE, rank,
                molpro::mpi::comm_global());
#endif
    auto range = distribution.range(molpro::mpi::rank_global());
    value_t value = 0;
    for (auto i = range.first; i < range.second; i++) {
      action_local_buffer[i - range.first] = 0;
      for (size_t j = 0; j < v.size(); j++)
        action_local_buffer[i - range.first] += matrix(i, j) * (m_vectors.front()[j] - 1);
      value += 0.5 * action_local_buffer[i - range.first] * (m_vectors.front()[i] - 1);
    }
#ifdef HAVE_MPI_H
    if (molpro::mpi::size_global() > 1)
      MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_DOUBLE, MPI_SUM, molpro::mpi::comm_global());
#endif
    return value;
  }

  void action(const CVecRef<container_t>& parameters, const VecRef<container_t>& actions) const override {
    for (size_t k = 0; k < parameters.size(); k++) {
      const auto& v = parameters[k].get();
      auto& a = actions[k].get();
      auto action_local_buffer_pointer = a.local_buffer();
      auto& action_local_buffer = *action_local_buffer_pointer;
      auto& distribution = a.distribution();
      assert(distribution.compatible(v.distribution()));
      assert(v.local_buffer().get()->data() <= m_vectors[k].data());
      assert(v.local_buffer().get()->data() + n > m_vectors[k].data());
#ifdef HAVE_MPI_H
      for (int rank = 0; rank < molpro::mpi::size_global(); rank++)
        MPI_Bcast((void*)(m_vectors[k].data() + distribution.range(rank).first),
                  distribution.range(rank).second - distribution.range(rank).first, MPI_DOUBLE, rank,
                  molpro::mpi::comm_global());
#endif
      auto range = distribution.range(molpro::mpi::rank_global());
      for (auto i = range.first; i < range.second; i++) {
        action_local_buffer[i - range.first] = 0;
        for (size_t j = 0; j < v.size(); j++)
          action_local_buffer[i - range.first] += matrix(i, j) * m_vectors[k][j];
      }
    }
  }
};
#endif // LINEARALGEBRA_EXAMPLES_EXAMPLEPROBLEM_H_
