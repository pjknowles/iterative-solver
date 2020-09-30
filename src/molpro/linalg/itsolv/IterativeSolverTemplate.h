#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_ITERATIVESOLVERTEMPLATE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_ITERATIVESOLVERTEMPLATE_H
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/Logger.h>
#include <molpro/linalg/itsolv/subspace/Matrix.h>
#include <molpro/linalg/itsolv/subspace/util.h>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace detail {

template <class R, class Q, class P>
void construct_solution(const std::vector<int>& working_set, std::vector<R>& params,
                        const std::vector<std::reference_wrapper<R>>& rparams,
                        const std::vector<std::reference_wrapper<Q>>& qparams,
                        const std::vector<std::reference_wrapper<P>>& pparams, size_t oR, size_t oQ, size_t oP,
                        const subspace::Matrix<double>& solutions, ArrayHandlers<R, Q, P>& handlers) {
  for (size_t i = 0; i < params.size(); ++i) {
    handlers.rr().fill(0, params.at(i));
  }
  for (size_t i = 0; i < working_set.size(); ++i) {
    handlers.rr().fill(0, params[i]);
    auto root = working_set[i];
    for (size_t j = 0; j < pparams.size(); ++j) {
      handlers.rp().axpy(solutions(root, oP + j), pparams.at(j), params.at(i));
    }
    for (size_t j = 0; j < qparams.size(); ++j) {
      handlers.rq().axpy(solutions(root, oQ + j), qparams.at(j), params.at(i));
    }
    for (size_t j = 0; j < rparams.size(); ++j) {
      handlers.rr().axpy(solutions(root, oR + j), rparams.at(j), params.at(i));
    }
  }
}
template <class R>
void normalise(const std::vector<int>& working_set, std::vector<R>& params, std::vector<R>& actions,
               array::ArrayHandler<R, R>& handler, Logger& logger) {
  for (size_t i = 0; i < working_set.size(); ++i) {
    auto dot = std::abs(handler.dot(params.at(i), params.at(i)));
    dot = std::sqrt(std::max(dot, decltype(dot)(0)));
    if (dot > 1.0e-14) {
      handler.scal(1. / dot, params.at(i));
      handler.scal(1. / dot, actions.at(i));
    } else {
      logger.msg("solution parameter's length is too small, dot = " + Logger::scientific(dot), Logger::Warn);
    }
  }
}

template <class R, typename T>
void construct_residual(const std::vector<int>& working_set, const std::vector<R>& solutions,
                        const std::vector<R>& actions, std::vector<std::reference_wrapper<R>>& residuals,
                        const std::vector<T>& eigvals, array::ArrayHandler<R, R>& handler) {
  assert(residuals.size() >= working_set.size());
  for (size_t i = 0; i < working_set.size(); ++i) {
    if (std::addressof(actions.at(i)) != std::addressof(residuals.at(i).get()))
      handler.copy(residuals.at(i), actions.at(i));
  }
  for (size_t i = 0; i < working_set.size(); ++i) {
    handler.axpy(-eigvals.at(working_set[i]), solutions.at(i), residuals.at(i));
  }
}

template <class R>
auto update_errors(const std::vector<int>& working_set, const std::vector<std::reference_wrapper<R>>& residual,
                   array::ArrayHandler<R, R>& handler) {
  auto errors = std::vector<double>(working_set.size());
  for (size_t i = 0; i < working_set.size(); ++i) {
    auto a = handler.dot(residual[i], residual[i]);
    errors[i] = std::sqrt(std::abs(a));
  }
  return errors;
}

} // namespace detail

/*!
 * @brief Implements common functionality of iterative solvers
 *
 * This is the trunk. It has a template of steps that all iterative solvers follow. Variations in implementation are
 * accepted as policies for managing the subspaces.
 *
 * Examples
 * ========
 * We are looking for n roots, but can only keep m < n roots in memory
 * @code{.cpp}
 * auto handlers = std::make_shared<ArrayHandlers<R, Q, P>>{};
 * auto solver = LinearEigensystemA<R, Q, P>{handlers};
 * solver.set_roots(n);
 * auto params = std::vector<R>{};
 * auto actions = std::vector<R>{};
 * initialize(params);
 * size_t n_work = 1; // number of working vectors
 * for (auto i = 0; i < max_it && n_work != 0; ++i){
 *   // calculate action of the matrix, one parameter at a time
 *   for (size_t i = 0; i < n_work; ++i){
 *     apply_matrix(params[i], actions[i]);
 *   }
 *   n_work = solver.add_vector(params, actions);
 *   solver.report();
 *   if (precondition_manually) {
 *     apply_preconditioner(params, actions, n_work);
 *   } else {
 *     solver.precondition(params[i], actions[i]);
 *   }
 *   n_work = solver.end_iteration(params, actions);
 * }
 * @endcode
 *
 *
 */
template <class Solver, class XS>
class IterativeSolverTemplate : public Solver {
public:
  using typename Solver::scalar_type;
  using typename Solver::value_type;
  using RS = typename XS::RS;
  using QS = typename XS::QS;
  using PS = typename XS::PS;
  using CS = typename XS::CS;
  using R = typename XS::R;
  using Q = typename XS::Q;
  using P = typename XS::P;
  template <typename T>
  using VecRef = std::vector<std::reference_wrapper<T>>;

  IterativeSolverTemplate() = delete;
  IterativeSolverTemplate(const IterativeSolverTemplate<Solver, XS>&) = delete;
  IterativeSolverTemplate(IterativeSolverTemplate<Solver, XS>&&) noexcept = default;
  IterativeSolverTemplate<Solver, XS>& operator=(const IterativeSolverTemplate<Solver, XS>&) = delete;
  IterativeSolverTemplate<Solver, XS>& operator=(IterativeSolverTemplate<Solver, XS>&&) noexcept = default;

  /*!
   * @brief Adds new parameters and corresponding action to the subspace and solves the corresponding problem.
   *
   * New algorithm
   * -------------
   *  - The parameters and action will be added to the Q space.
   *  - The P space, I'm not sure yet.
   *  - The R space just wraps new parameters.
   *  - The Q space is updated using R space (e.g. for linear eigenvalue problem just copies params and clears R space).
   *  - The S space contains best solutions so far (empty on the start).
   *  - The full X subspace is formed.
   *  - The X space is checked for conditioning and vectors may be removed depending on implementation.
   *  - The subspace problem is solved and solutions are stored in the S space.
   *  - The working set of vectors is made up of vectors with largest errors that are not converged.
   *
   * @param parameters new parameters for the R space
   * @param action corresponding action
   * @param parametersP ???
   * @return
   */
  size_t add_vector(std::vector<R>& parameters, std::vector<R>& action, std::vector<P>& parametersP) override {
    using subspace::util::wrap;
    assert(parameters.size() >= m_working_set.size());
    assert(action.size() >= m_working_set.size());
    m_logger->msg("IterativeSolverTemplate::add_vector  iteration = " + std::to_string(m_stats->iterations),
                  Logger::Trace);
    m_logger->msg("IterativeSolverTemplate::add_vector  size of {params, actions, working_set} = " +
                      std::to_string(parameters.size()) + ", " + std::to_string(action.size()) + ", " +
                      std::to_string(m_working_set.size()) + ", ",
                  Logger::Debug);
    m_rspace.update(parameters, action, *static_cast<Solver*>(this));
    m_qspace.update(m_rspace, m_cspace, *static_cast<Solver*>(this));
    m_xspace.build_subspace(m_rspace, m_qspace, m_pspace, m_cspace);
    m_xspace.check_conditioning(m_rspace, m_qspace, m_pspace);
    m_xspace.solve(*static_cast<Solver*>(this));
    auto& dummy = m_rspace.dummy(parameters.size());
    auto wdummy = wrap(dummy);
    // FIXME Constructing solution, calculating errors and updating the working space all have to be done in one go
    detail::construct_solution(m_working_set, parameters, wrap(m_rspace.params()), m_qspace.params(), m_pspace.params(),
                               m_xspace.dimensions().oR, m_xspace.dimensions().oQ, m_xspace.dimensions().oP,
                               m_xspace.solutions(), *m_handlers);
    detail::construct_solution(m_working_set, action, wrap(m_rspace.actions()), m_qspace.actions(), m_pspace.actions(),
                               m_xspace.dimensions().oR, m_xspace.dimensions().oQ, m_xspace.dimensions().oP,
                               m_xspace.solutions(), *m_handlers);
    detail::normalise(m_working_set, parameters, action, m_handlers->rr(), *m_logger);
    detail::construct_residual(m_working_set, parameters, action, wdummy, m_xspace.eigenvalues(), m_handlers->rr());
    m_errors = detail::update_errors(m_working_set, wdummy, m_handlers->rr());
    m_logger->msg("add_vector::errors = ", begin(m_errors), end(m_errors), Logger::Trace);
    update_working_set();
    for (size_t i = 0; i < action.size(); ++i) {
      m_handlers->rr().copy(action.at(i), dummy.at(i));
    }
    m_stats->iterations++;
    return m_working_set.size();
  };

  size_t add_vector(std::vector<R>& parameters, std::vector<R>& action) override {
    auto p = std::vector<P>{};
    return add_vector(parameters, action, p);
  }

  size_t add_p(std::vector<P>& Pvectors, const value_type* PP, std::vector<R>& parameters, std::vector<R>& action,
               std::vector<P>& parametersP) override {
    return 0;
  };

  /*!
   * @brief Copy the solutions from the S space
   *
   * @param roots
   * @param parameters
   * @param residual
   */
  void solution(const std::vector<int>& roots, std::vector<R>& parameters, std::vector<R>& residual) override{};

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
  void set_n_roots(size_t roots) override { m_nroots = roots; }
  const std::vector<scalar_type>& errors() const override { return m_errors; }
  const Statistics& statistics() const override { return *m_stats; }

  void report() const override {
    molpro::cout << "iteration " << m_stats->iterations;
    if (not m_errors.empty()) {
      auto it_max_error = std::max_element(m_errors.cbegin(), m_errors.cend());
      if (n_roots() > 1)
        molpro::cout << ", error[" << std::distance(m_errors.cbegin(), it_max_error) << "] = ";
      else
        molpro::cout << ", error = ";
      molpro::cout << *it_max_error << std::endl;
    }
  }

protected:
  IterativeSolverTemplate(RS rspace, QS qspace, PS pspace, CS sspace, XS xspace,
                          std::shared_ptr<ArrayHandlers<R, Q, P>> handlers, std::shared_ptr<Statistics> stats,
                          std::shared_ptr<Logger> logger)
      : m_handlers(std::move(handlers)), m_rspace(std::move(rspace)), m_qspace(std::move(qspace)),
        m_pspace(std::move(pspace)), m_cspace(std::move(sspace)), m_xspace(std::move(xspace)),
        m_stats(std::move(stats)), m_logger(std::move(logger)) {}

  //! Updates working sets and adds current solutions to the S space
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
    m_logger->msg("update_working_set::converged_roots = ", begin(converged_roots), end(converged_roots),
                  Logger::Debug);
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
  CS m_cspace;
  XS m_xspace;
  std::vector<double> m_errors;
  std::vector<int> m_working_set;
  size_t m_nroots{0};
  double m_convergence_threshold{1.0e-10}; //!< errors less than this mark a converged solution
  std::shared_ptr<Statistics> m_stats;
  std::shared_ptr<Logger> m_logger;
};

} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_ITERATIVESOLVERTEMPLATE_H
