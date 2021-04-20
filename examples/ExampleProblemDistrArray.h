#ifndef LINEARALGEBRA_EXAMPLES_EXAMPLEPROBLEM_H_
#define LINEARALGEBRA_EXAMPLES_EXAMPLEPROBLEM_H_
#include <molpro/linalg/array/DistrArraySpan.h>
#include <molpro/linalg/array/util/Distribution.h>
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/mpi.h>
#include <vector>

class ExampleProblem : public molpro::linalg::itsolv::Problem<molpro::linalg::array::DistrArraySpan> {

protected:
  double matrix(int i, int j) const { return i == j ? i + 1.5 : 0.001 * ((i + j) % n); }

public:
  using Pvector = std::map<size_t, double>;
  using Problem::container_t;
  using Problem::value_t;
  const size_t n;
  std::vector<std::vector<double>> m_vectors;
  ExampleProblem(int n = 10, size_t nvector = 1) : n(n) {
    for (size_t i = 0; i < nvector; i++) {
      m_vectors.emplace_back(n);
    }
  }

  bool diagonals(container_t& d) const override {
    auto range = d.distribution().range(molpro::mpi::rank_global());
    auto d_local_buffer = d.local_buffer();
    for (auto i = range.first; i < range.second; i++)
      (*d_local_buffer)[i - range.first] = matrix(i, i);
    return true;
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

  std::vector<double> pp_action_matrix(const std::vector<Pvector>& pparams) const override {
    std::vector<double> result(pparams.size() * pparams.size(), 0);
    size_t ij = 0;
    for (const auto& pi : pparams)
      for (const auto& pj : pparams) {
        for (const auto& pie : pi)
          for (const auto& pje : pj)
            result[ij] = matrix(pje.first, pie.first) * pje.second * pie.second;
        ij++;
      }
    return result;
  }

  void p_action(const std::vector<std::vector<value_t>>& p_coefficients, const CVecRef<Pvector>& pparams,
                const VecRef<container_t>& actions) const override {
    for (size_t k = 0; k < p_coefficients.size(); k++) {
      auto& a = actions[k].get();
      auto action_local_buffer_pointer = a.local_buffer();
      auto& action_local_buffer = *action_local_buffer_pointer;
      auto range = a.distribution().range(molpro::mpi::rank_global());
      for (size_t pindex = 0; pindex < pparams.size(); pindex++) {
        const auto& pi = pparams[pindex].get();
        for (const auto& pie : pi) {
          auto coeff = pie.second * p_coefficients[k][pindex];
          for (auto i = range.first; i < range.second; i++)
            action_local_buffer[i - range.first] += matrix(i, pie.first) * coeff;
        }
      }
    }
  }
};
#endif // LINEARALGEBRA_EXAMPLES_EXAMPLEPROBLEM_H_
