#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_ITERATIVESOLVERTEMPLATE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_ITERATIVESOLVERTEMPLATE_H
#include <molpro/linalg/itsolv/IterativeSolver.h>

namespace molpro {
namespace linalg {
namespace itsolv {

namespace detail {} // namespace detail

/*!
 * @brief Implements common functionality of iterative solvers
 *
 * This is the trunk. It has a template of steps that all iterative solvers follow. Variations in implementation are
 * accepted as policies for managing the subspaces.
 *
 */
template <class Solver, class RSpace, class QSpace, class PSpace, class XSpace>
class IterativeSolverTemplate : public Solver {
public:
  using typename Solver::scalar_type;
  using R = typename RSpace::R;
  using Q = typename QSpace::Q;
  using P = typename PSpace::P;

  void add_vector(std::vector<R>& parameters, std::vector<R>& action, std::vector<P>& parametersP) override {
    m_rspace.update(parameters, action, *this); // allows for passage of extra data from the particular Solver
    m_qspace.add(m_rspace, *this);
    m_xspace.build_subspace(m_rspace, m_qspace, m_pspace);
    m_xspace.check_conditioning(m_rspace, m_qspace, m_pspace);
    m_xspace.solve();
    construct_solution(parameters);
    construct_residual(parameters, action);
    update_errors(action);
    construct_action(action);
    update_working_set();
    construct_residual(parameters, action);
  };

  void solution(const std::vector<int>& roots, std::vector<R>& parameters, std::vector<R>& residual) override {
    auto working_set_save = m_working_set;
    m_working_set = roots;
    m_xspace.build_subspace(m_rspace, m_qspace, m_pspace);
    m_xspace.solve();
    construct_solution(parameters);
    construct_residual(residual);
    m_working_set = working_set_save;
  };

  void solution(const std::vector<int>& roots, std::vector<R>& parameters, std::vector<R>& residual,
                std::vector<P>& parametersP) override {}

  std::vector<size_t> suggest_p(const std::vector<R>& solution, const std::vector<R>& residual, size_t maximumNumber,
                                double threshold) override {}

  const std::vector<int>& working_set() const override { return m_working_set; }
  const std::vector<scalar_type>& errors() const override { return m_errors; }
  const Statistics& statistics() const override { return *m_stats; }

protected:
  //! Recalculates the error vector
  void update_errors(std::vector<R>& residual) {}

  //! Updates working sets and adds any converged solution to the q space
  void update_working_set() {}

  void construct_solution(std::vector<typename RSpace::R>& solution) {}

  void construct_residual(const std::vector<typename RSpace::R>& solution, std::vector<typename RSpace::R>& residual) {}

  void construct_action(std::vector<typename RSpace::R>& action) {}

  std::shared_ptr<ArrayHandlers<typename RSpace::R, typename QSpace::Q, typename PSpace::P>> m_handlers;
  RSpace m_rspace;
  QSpace m_qspace;
  PSpace m_pspace;
  XSpace m_xspace;
  std::vector<double> m_errors;
  std::vector<int> m_working_set;
  std::shared_ptr<Statistics> m_stats;
};

} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_ITERATIVESOLVERTEMPLATE_H
