#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_ITERATIVESOLVERTEMPLATE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_ITERATIVESOLVERTEMPLATE_H
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/Logger.h>
#include <molpro/linalg/itsolv/subspace/Matrix.h>
#include <molpro/linalg/itsolv/subspace/util.h>
#include <molpro/linalg/itsolv/wrap.h>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace detail {

std::vector<std::pair<size_t, size_t>> parameter_batches(const size_t nsol, const size_t nparam) {
  auto n_batch = nsol / nparam + (nsol == nparam ? 0 : 1);
  auto batches = std::vector<std::pair<size_t, size_t>>{};
  for (size_t ib = 0, start_sol = 0, end_sol = 0; ib < n_batch; ++ib, start_sol = end_sol) {
    end_sol = (start_sol + nparam > nsol ? nsol - start_sol : start_sol + nparam);
    batches.emplace_back(start_sol, end_sol);
  }
  return batches;
}

template <class R, class Q, class P>
void construct_solution(std::vector<R>& params, const std::vector<unsigned int>& roots,
                        const subspace::Matrix<double>& solutions,
                        const std::vector<std::reference_wrapper<P>>& pparams,
                        const std::vector<std::reference_wrapper<Q>>& qparams,
                        const std::vector<std::reference_wrapper<Q>>& dparams, size_t oP, size_t oQ, size_t oD,
                        ArrayHandlers<R, Q, P>& handlers) {
  assert(params.size() >= roots.size());
  for (size_t i = 0; i < roots.size(); ++i) {
    handlers.rr().fill(0, params.at(i));
  }
  for (size_t i = 0; i < roots.size(); ++i) {
    auto root = roots[i];
    for (size_t j = 0; j < pparams.size(); ++j) {
      handlers.rp().axpy(solutions(root, oP + j), pparams.at(j), params.at(i));
    }
    for (size_t j = 0; j < qparams.size(); ++j) {
      handlers.rq().axpy(solutions(root, oQ + j), qparams.at(j), params.at(i));
    }
    for (size_t j = 0; j < dparams.size(); ++j) {
      handlers.rq().axpy(solutions(root, oD + j), dparams.at(j), params.at(i));
    }
  }
}

template <class R, class P>
void remove_p_component(std::vector<R>& params, const std::vector<unsigned int>& roots,
                        const subspace::Matrix<double>& solutions,
                        const std::vector<std::reference_wrapper<P>>& pparams, size_t oP,
                        array::ArrayHandler<R, P>& handler) {
  assert(params.size() >= roots.size());
  for (size_t i = 0; i < roots.size(); ++i) {
    auto root = roots[i];
    for (size_t j = 0; j < pparams.size(); ++j) {
      handler.axpy(-solutions(root, oP + j), pparams.at(j), params.at(i));
    }
  }
}

template <class P, typename T>
std::vector<P> construct_pspace_vector(const std::vector<unsigned int>& roots, const subspace::Matrix<T>& solutions,
                                       const std::vector<std::reference_wrapper<P>>& pparams, const size_t oP,
                                       array::ArrayHandler<P, P>& handler) {
  auto pvecs = std::vector<P>{};
  if (!pparams.empty()) {
    for (size_t i = 0; i < roots.size(); ++i) {
      auto root = roots[i];
      pvecs.emplace_back(handler.copy(pparams.front()));
      handler.scal(solutions(root, oP), pvecs.back());
      for (size_t j = 1; j < pparams.size(); ++j) {
        handler.axpy(solutions(root, oP + j), pparams[j], pvecs.back());
      }
    }
  }
  return pvecs;
}

template <class R>
void normalise(const size_t n_roots, std::vector<R>& params, std::vector<R>& actions,
               array::ArrayHandler<R, R>& handler, Logger& logger) {
  assert(params.size() >= n_roots && actions.size() >= n_roots);
  for (size_t i = 0; i < n_roots; ++i) {
    auto dot = handler.dot(params.at(i), params.at(i));
    dot = std::sqrt(std::abs(dot));
    if (dot > 1.0e-14) {
      handler.scal(1. / dot, params.at(i));
      handler.scal(1. / dot, actions.at(i));
    } else {
      logger.msg("solution parameter's length is too small, dot = " + Logger::scientific(dot), Logger::Warn);
    }
  }
}

template <class R, typename T>
void construct_residual(const std::vector<unsigned int>& roots, const std::vector<T>& eigvals,
                        const std::vector<R>& params, std::vector<R>& actions, array::ArrayHandler<R, R>& handler) {
  assert(params.size() >= roots.size());
  for (size_t i = 0; i < roots.size(); ++i) {
    handler.axpy(-eigvals.at(roots[i]), params.at(i), actions.at(i));
  }
}

template <class R, typename T>
void update_errors(std::vector<T>& errors, const std::vector<R>& residual, array::ArrayHandler<R, R>& handler) {
  assert(residual.size() >= errors.size());
  for (size_t i = 0; i < errors.size(); ++i) {
    auto a = handler.dot(residual[i], residual[i]);
    errors[i] = std::sqrt(std::abs(a));
  }
}

template <typename T>
std::vector<unsigned int> select_working_set(const size_t nw, const std::vector<T>& errors, const T threshold) {
  auto ordered_errors = std::multimap<T, size_t, std::greater<T>>{};
  for (size_t i = 0; i < errors.size(); ++i) {
    if (errors[i] > threshold)
      ordered_errors.emplace(errors[i], i);
  }
  auto working_set = std::vector<unsigned int>{};
  auto end = (ordered_errors.size() < nw ? ordered_errors.end() : next(begin(ordered_errors), nw));
  std::transform(begin(ordered_errors), end, std::back_inserter(working_set), [](const auto& el) { return el.second; });
  return working_set;
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
template <class Solver, class XS, class SubspaceSolver>
class IterativeSolverTemplate : public Solver {
public:
  using typename Solver::fapply_on_p_type;
  using typename Solver::scalar_type;
  using typename Solver::value_type;
  using R = typename XS::R;
  using Q = typename XS::Q;
  using P = typename XS::P;

  IterativeSolverTemplate() = delete;
  IterativeSolverTemplate(const IterativeSolverTemplate<Solver, XS, SubspaceSolver>&) = delete;
  IterativeSolverTemplate(IterativeSolverTemplate<Solver, XS, SubspaceSolver>&&) noexcept = default;
  IterativeSolverTemplate<Solver, XS, SubspaceSolver>&
  operator=(const IterativeSolverTemplate<Solver, XS, SubspaceSolver>&) = delete;
  IterativeSolverTemplate<Solver, XS, SubspaceSolver>&
  operator=(IterativeSolverTemplate<Solver, XS, SubspaceSolver>&&) noexcept = default;

protected:
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
   * @param pparams P space components of the working set solutions
   * @return
   */
  size_t add_vector(std::vector<R>& parameters, std::vector<R>& action, std::vector<P>& pparams,
                    fapply_on_p_type& apply_p) {
    assert(parameters.size() == action.size());
    m_logger->msg("IterativeSolverTemplate::add_vector  iteration = " + std::to_string(m_stats->iterations) +
                      ", apply_p = " + std::to_string(bool(apply_p)),
                  Logger::Trace);
    m_logger->msg("IterativeSolverTemplate::add_vector  size of {params, actions, working_set} = " +
                      std::to_string(parameters.size()) + ", " + std::to_string(action.size()) + ", " +
                      std::to_string(m_working_set.size()) + ", ",
                  Logger::Debug);
    auto nW = std::min(m_working_set.size(), parameters.size());
    auto wparams = cwrap<R>(begin(parameters), begin(parameters) + nW);
    auto wactions = cwrap<R>(begin(action), begin(action) + nW);
    // FIXME finish constructing D space action
    m_xspace.update_qspace(wparams, wactions);
    m_subspace_solver.solve(m_xspace, n_roots());
    auto nsol = m_subspace_solver.size();
    std::vector<std::pair<Q, Q>> temp_solutions{};
    for (const auto& batch : detail::parameter_batches(nsol, parameters.size())) {
      size_t start_sol, end_sol;
      std::tie(start_sol, end_sol) = batch;
      auto roots = std::vector<unsigned int>(end_sol - start_sol);
      std::iota(begin(roots), end(roots), start_sol);
      detail::construct_solution(parameters, roots, m_subspace_solver.solutions(), m_xspace.paramsp(),
                                 m_xspace.paramsq(), m_xspace.paramsd(), m_xspace.dimensions().oP,
                                 m_xspace.dimensions().oQ, m_xspace.dimensions().oD, *m_handlers);
      detail::construct_solution(action, roots, m_subspace_solver.solutions(), m_xspace.actionsp(), m_xspace.actionsq(),
                                 m_xspace.actionsd(), m_xspace.dimensions().oP, m_xspace.dimensions().oQ,
                                 m_xspace.dimensions().oD, *m_handlers);
      auto pvectors = detail::construct_pspace_vector(roots, m_subspace_solver.solutions(), m_xspace.paramsp(),
                                                      m_xspace.dimensions().oP, m_handlers->pp());
      detail::normalise(roots.size(), parameters, action, m_handlers->rr(), *m_logger);
      if (apply_p) {
        auto waction = wrap(action);
        apply_p(pvectors, waction);
      } else {
        detail::remove_p_component(parameters, roots, m_subspace_solver.solutions(), m_xspace.paramsp(),
                                   m_xspace.dimensions().oP, m_handlers->rp());
        detail::remove_p_component(action, roots, m_subspace_solver.solutions(), m_xspace.actionsp(),
                                   m_xspace.dimensions().oP, m_handlers->rp());
      }
      detail::construct_residual(roots, m_subspace_solver.eigenvalues(), parameters, action, m_handlers->rr());
      auto errors = std::vector<scalar_type>(roots.size(), 0);
      detail::update_errors(errors, action, m_handlers->rr());
      for (size_t i = 0; i < roots.size(); ++i)
        temp_solutions.emplace_back(m_handlers->qr().copy(parameters[i]), m_handlers->qr().copy(action[i]));
      m_subspace_solver.set_error(roots, errors);
    }
    m_errors = m_subspace_solver.errors();
    m_working_set = detail::select_working_set(parameters.size(), m_errors, m_convergence_threshold);
    for (size_t i = 0; i < m_working_set.size(); ++i) {
      auto root = m_working_set[i];
      m_handlers->rq().copy(parameters[i], temp_solutions.at(root).first);
      m_handlers->rq().copy(action[i], temp_solutions.at(root).second);
    }
    pparams = detail::construct_pspace_vector(m_working_set, m_subspace_solver.solutions(), m_xspace.paramsp(),
                                              m_xspace.dimensions().oP, m_handlers->pp());
    m_logger->msg("add_vector::errors = ", begin(m_errors), end(m_errors), Logger::Trace);
    return m_working_set.size();
  }

public:
  size_t add_vector(std::vector<R>& parameters, std::vector<R>& action, fapply_on_p_type& apply_p) override {
    auto pparams = std::vector<P>{};
    return add_vector(parameters, action, pparams, apply_p);
  }

  size_t add_vector(std::vector<R>& parameters, std::vector<R>& action, std::vector<P>& pparams) override {
    auto apply_p = fapply_on_p_type{};
    return add_vector(parameters, action, pparams, apply_p);
  }

  size_t add_vector(std::vector<R>& parameters, std::vector<R>& action) override {
    auto pparams = std::vector<P>{};
    auto apply_p = fapply_on_p_type{};
    return add_vector(parameters, action, pparams, apply_p);
  }

  // TODO Add P space and solve the subspace problem. This is same as add_vector, but without updating the Q space
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
  void solution(const std::vector<unsigned int>& roots, std::vector<R>& parameters,
                std::vector<R>& residual) override{};

  void solution(const std::vector<unsigned int>& roots, std::vector<R>& parameters, std::vector<R>& residual,
                std::vector<P>& parametersP) override {
    solution(roots, parameters, residual);
  }

  std::vector<size_t> suggest_p(const std::vector<R>& solution, const std::vector<R>& residual, size_t maximumNumber,
                                double threshold) override {
    return {};
  }

  const std::vector<unsigned int>& working_set() const override { return m_working_set; }

  size_t n_roots() const override { return m_nroots; }

  void set_n_roots(size_t roots) override {
    m_nroots = roots;
    m_working_set.resize(roots);
    std::iota(begin(m_working_set), end(m_working_set), (unsigned int)0);
  }

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
  IterativeSolverTemplate(XS xspace, SubspaceSolver solver, std::shared_ptr<ArrayHandlers<R, Q, P>> handlers,
                          std::shared_ptr<Statistics> stats, std::shared_ptr<Logger> logger)
      : m_handlers(std::move(handlers)), m_xspace(std::move(xspace)), m_subspace_solver(std::move(solver)),
        m_stats(std::move(stats)), m_logger(std::move(logger)) {}

  std::shared_ptr<ArrayHandlers<R, Q, P>> m_handlers;
  XS m_xspace;
  SubspaceSolver m_subspace_solver;
  std::vector<double> m_errors;
  std::vector<unsigned int> m_working_set;
  size_t m_nroots{0};
  double m_convergence_threshold{1.0e-10}; //!< errors less than this mark a converged solution
  std::shared_ptr<Statistics> m_stats;
  std::shared_ptr<Logger> m_logger;
};

} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_ITERATIVESOLVERTEMPLATE_H
