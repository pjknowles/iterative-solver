#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_ITERATIVESOLVERTEMPLATE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_ITERATIVESOLVERTEMPLATE_H
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/subspace/Matrix.h>
#include <molpro/linalg/itsolv/subspace/util.h>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace detail {

template <class R, class Q, class P>
void construct_solution(const std::vector<int>& working_set, std::vector<std::reference_wrapper<R>>& params,
                        const std::vector<std::reference_wrapper<R>>& dummy,
                        const std::vector<std::reference_wrapper<Q>>& qparams,
                        const std::vector<std::reference_wrapper<P>>& pparams, size_t oR, size_t oQ, size_t oP,
                        const subspace::Matrix<double>& solutions, ArrayHandlers<R, Q, P>& handlers) {
  assert(dummy.size() >= working_set.size());
  for (size_t i = 0; i < working_set.size(); ++i) {
    handlers.rr().copy(dummy.at(i), params.at(i));
  }
  for (size_t i = 0; i < working_set.size(); ++i) {
    handlers.rr().fill(0, params[i]);
    auto solution = solutions.col(working_set[i]);
    for (size_t j = 0; j < pparams.size(); ++j) {
      handlers.rp().axpy(solution(oP + j, 0), pparams.at(j), params.at(i));
    }
    for (size_t j = 0; j < qparams.size(); ++j) {
      handlers.rq().axpy(solution(oQ + j, 0), qparams.at(j), params.at(i));
    }
    for (size_t j = 0; j < dummy.size(); ++j) {
      handlers.rr().axpy(solution(oR + j, 0), dummy.at(j), params.at(i));
    }
  }
}

template <class R, typename T>
void construct_residual(const std::vector<int>& working_set, const std::vector<std::reference_wrapper<R>>& solutions,
                        const std::vector<std::reference_wrapper<R>>& actions,
                        std::vector<std::reference_wrapper<R>>& residuals, const std::vector<T>& eigvals,
                        array::ArrayHandler<R, R>& handler) {
  assert(residuals.size() >= working_set.size());
  for (size_t i = 0; i < working_set.size(); ++i) {
    handler.rr().copy(residuals.at(i), actions.at(i));
  }
  for (size_t i = 0; i < working_set.size(); ++i) {
    handler.axpy(-eigvals.at(working_set[i]), solutions.at(i), residuals.at(i));
  }
}

template <class R>
auto update_errors(const std::vector<std::reference_wrapper<R>>& residual, array::ArrayHandler<R, R>& handler) {
  auto errors = std::vector<double>(residual.size());
  for (size_t i = 0; i < residual.size(); ++i) {
    auto a = handler->rr().dot(residual[i], residual[i]);
    errors[i] = std::sqrt(std::abs(a));
  }
}

} // namespace detail

/*!
 * @brief Implements common functionality of iterative solvers
 *
 * This is the trunk. It has a template of steps that all iterative solvers follow. Variations in implementation are
 * accepted as policies for managing the subspaces.
 *
 */
template <class Solver, class XS>
class IterativeSolverTemplate : public Solver {
public:
  using typename Solver::scalar_type;
  using RS = typename XS::RS;
  using QS = typename XS::QS;
  using PS = typename XS::PS;
  using R = typename Solver::R;
  using Q = typename Solver::Q;
  using P = typename Solver::P;
  template <typename T>
  using VecRef = std::vector<std::reference_wrapper<T>>;

  IterativeSolverTemplate() = delete;
  IterativeSolverTemplate(const IterativeSolverTemplate<Solver, XS>&) = delete;
  IterativeSolverTemplate(IterativeSolverTemplate<Solver, XS>&&) noexcept = default;
  IterativeSolverTemplate<Solver, XS>& operator=(const IterativeSolverTemplate<Solver, XS>&) = delete;
  IterativeSolverTemplate<Solver, XS>& operator=(IterativeSolverTemplate<Solver, XS>&&) noexcept = default;

  void add_vector(std::vector<R>& parameters, std::vector<R>& action, std::vector<P>& parametersP) override {
    using subspace::util::wrap;
    m_rspace.update(parameters, action, *static_cast<Solver*>(this));
    m_working_set.clear();
    std::copy(begin(m_rspace.working_set), end(m_rspace.working_set()), std::back_inserter(m_working_set));
    if (m_nroots == 0)
      m_nroots = m_working_set.size();
    m_qspace.update(m_rspace, *static_cast<Solver*>(this));
    m_xspace.build_subspace(m_rspace, m_qspace, m_pspace);
    m_xspace.check_conditioning(m_rspace, m_qspace, m_pspace);
    m_xspace.solve(*static_cast<Solver*>(this));
    auto& dummy = m_rspace.dummy(m_working_set.size());
    detail::construct_solution(m_working_set, m_rspace.params(), wrap(dummy), m_qspace.params(), m_pspace.params(),
                               m_xspace.dimensions().oR, m_xspace.dimensions().oQ, m_xspace.dimensions().oP,
                               m_xspace.solutions(), *m_handlers);
    detail::construct_solution(m_working_set, m_rspace.actions(), wrap(dummy), m_qspace.actions(), m_pspace.actions(),
                               m_xspace.dimensions().oR, m_xspace.dimensions().oQ, m_xspace.dimensions().oP,
                               m_xspace.solutions(), *m_handlers);
    detail::construct_residual(m_working_set, m_rspace.params(), m_rspace.actions(), wrap(dummy),
                               m_xspace.eigenvalues(), m_handlers->rr());
    m_errors = detail::update_errors(wrap(dummy), m_handlers->rr());
    update_working_set();
    for (size_t i = 0; i < m_working_set.size(); ++i)
      m_handlers->rr().copy(m_rspace.actions().at(i), dummy.at(i));
    dummy.clear();
  };

  // FIXME I don't fully understand what this is supposed to be doing and what the input is
  void solution(const std::vector<int>& roots, std::vector<R>& parameters, std::vector<R>& residual) override {
    auto working_set_save = m_working_set;
    m_working_set = roots;
    m_xspace.build_subspace(m_rspace, m_qspace, m_pspace);
    m_xspace.solve(*static_cast<Solver*>(this));
    //    construct_solution(parameters);
    //    construct_residual(residual);
    m_working_set = working_set_save;
  };

  void solution(const std::vector<int>& roots, std::vector<R>& parameters, std::vector<R>& residual,
                std::vector<P>& parametersP) override {
    solution(roots, parameters, residual);
  }

  std::vector<size_t> suggest_p(const std::vector<R>& solution, const std::vector<R>& residual, size_t maximumNumber,
                                double threshold) override {
    return {};
  }

  const std::vector<int>& working_set() const override { return m_working_set; }
  size_t n_roots() const override { return m_nroots; }
  const std::vector<scalar_type>& errors() const override { return m_errors; }
  const Statistics& statistics() const override { return *m_stats; }

protected:
  IterativeSolverTemplate(std::shared_ptr<RS> rspace, std::shared_ptr<QS> qspace, std::shared_ptr<PS> pspace,
                          std::shared_ptr<XS> xspace, std::shared_ptr<ArrayHandlers<R, Q, P>> handlers,
                          std::shared_ptr<Statistics> stats)
      : m_handlers(std::move(handlers)), m_rspace(std::move(rspace)), m_qspace(std::move(qspace)),
        m_pspace(std::move(pspace)), m_xspace(std::move(xspace)), m_stats(std::move(stats)) {}

  //! Updates working sets and adds any converged solution to the q space
  void update_working_set() {
    auto ind_still_a_working_param = std::vector<size_t>{};
    auto converged_roots = std::vector<size_t>{};
    auto converged_params = VecRef<R>{};
    auto converged_actions = VecRef<R>{};
    for (size_t i = 0; i < m_working_set.size(); ++i) {
      if (m_errors.at(i) < m_convergence_threshold) {
        converged_roots.emplace_back(m_working_set[i]);
        converged_params.emplace_back(m_rspace.params().at(i));
        converged_actions.emplace_back(m_rspace.actions().at(i));
      } else {
        ind_still_a_working_param.emplace_back(i);
      }
    }
    m_qspace.add_converged(converged_params, converged_actions, converged_roots);
    m_rspace.update_working_set(ind_still_a_working_param);
    auto& new_working_set = m_rspace.working_set();
    m_working_set.resize(new_working_set.size());
    std::copy(begin(new_working_set), end(new_working_set), begin(m_working_set));
  }

  std::shared_ptr<ArrayHandlers<R, Q, P>> m_handlers;
  RS m_rspace;
  QS m_qspace;
  PS m_pspace;
  XS m_xspace;
  std::vector<double> m_errors;
  std::vector<int> m_working_set;
  size_t m_nroots{0};
  double m_convergence_threshold{1.0e-10}; //!< errors less than this mark a converged solution
  std::shared_ptr<Statistics> m_stats;
};

} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_ITERATIVESOLVERTEMPLATE_H
