#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_ITERATIVESOLVERTEMPLATE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_ITERATIVESOLVERTEMPLATE_H
#include <iostream>
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/Logger.h>
#include <molpro/linalg/itsolv/subspace/Matrix.h>
#include <molpro/linalg/itsolv/subspace/SubspaceSolverI.h>
#include <molpro/linalg/itsolv/subspace/XSpaceI.h>
#include <molpro/linalg/itsolv/subspace/util.h>
#include <molpro/linalg/itsolv/wrap.h>

namespace molpro::linalg::itsolv {
namespace detail {

std::vector<std::pair<size_t, size_t>> parameter_batches(const size_t nsol, const size_t nparam) {
  auto n_batch = nsol / nparam + (nsol % nparam ? 1 : 0);
  auto batches = std::vector<std::pair<size_t, size_t>>{};
  for (size_t ib = 0, start_sol = 0, end_sol = 0; ib < n_batch; ++ib, start_sol = end_sol) {
    end_sol = std::min(start_sol + nparam, nsol);
    batches.emplace_back(start_sol, end_sol);
  }
  return batches;
}

template <class R, class Q, class P>
void construct_solution(const VecRef<R>& params, const std::vector<int>& roots,
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
void remove_p_component(const VecRef<R>& params, const std::vector<int>& roots,
                        const subspace::Matrix<double>& solutions, const CVecRef<P>& pparams, size_t oP,
                        array::ArrayHandler<R, P>& handler) {
  assert(params.size() >= roots.size());
  for (size_t i = 0; i < roots.size(); ++i) {
    auto root = roots[i];
    for (size_t j = 0; j < pparams.size(); ++j) {
      handler.axpy(-solutions(root, oP + j), pparams.at(j), params.at(i));
    }
  }
}

template <typename T>
std::vector<std::vector<T>> construct_vectorP(const std::vector<int>& roots, const subspace::Matrix<T>& solutions,
                                              const size_t oP, const size_t nP) {
  auto vectorP = std::vector<std::vector<T>>{};
  for (auto root : roots) {
    vectorP.emplace_back();
    for (size_t j = 0; j < nP; ++j)
      vectorP.back().push_back(solutions(root, oP + j));
  }
  return vectorP;
}

template <class R>
void normalise(const size_t n_roots, const VecRef<R>& params, const VecRef<R>& actions,
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
void construct_residual(const std::vector<int>& roots, const std::vector<T>& eigvals, const CVecRef<R>& params,
                        const VecRef<R>& actions, array::ArrayHandler<R, R>& handler) {
  assert(params.size() >= roots.size());
  for (size_t i = 0; i < roots.size(); ++i) {
    handler.axpy(-eigvals.at(roots[i]), params.at(i), actions.at(i));
  }
}

template <class R, typename T>
void update_errors(std::vector<T>& errors, const CVecRef<R>& residual, array::ArrayHandler<R, R>& handler) {
  assert(residual.size() >= errors.size());
  for (size_t i = 0; i < errors.size(); ++i) {
    auto a = handler.dot(residual[i], residual[i]);
    errors[i] = std::sqrt(std::abs(a));
  }
}

template <typename T>
std::vector<int> select_working_set(const size_t nw, const std::vector<T>& errors, const T threshold) {
  auto ordered_errors = std::multimap<T, size_t, std::greater<T>>{};
  for (size_t i = 0; i < errors.size(); ++i) {
    if (errors[i] > threshold)
      ordered_errors.emplace(errors[i], i);
  }
  auto working_set = std::vector<int>{};
  auto end = (ordered_errors.size() < nw ? ordered_errors.end() : next(begin(ordered_errors), nw));
  std::transform(begin(ordered_errors), end, std::back_inserter(working_set), [](const auto& el) { return el.second; });
  return working_set;
}

} // namespace detail

/*!
 * @brief Implements functionality common to all iterative solvers
 */
template <template <class, class, class> class Solver, class R, class Q, class P>
class IterativeSolverTemplate : public Solver<R, Q, P> {
public:
  using typename Solver<R, Q, P>::fapply_on_p_type;
  using typename Solver<R, Q, P>::scalar_type;
  using typename Solver<R, Q, P>::value_type;
  using typename Solver<R, Q, P>::VectorP;

  IterativeSolverTemplate() = delete;
  IterativeSolverTemplate(const IterativeSolverTemplate<Solver, R, Q, P>&) = delete;
  IterativeSolverTemplate(IterativeSolverTemplate<Solver, R, Q, P>&&) noexcept = default;
  IterativeSolverTemplate<Solver, R, Q, P>& operator=(const IterativeSolverTemplate<Solver, R, Q, P>&) = delete;
  IterativeSolverTemplate<Solver, R, Q, P>& operator=(IterativeSolverTemplate<Solver, R, Q, P>&&) noexcept = default;

protected:
  /*!
   * @brief Adds new parameters and corresponding action to the subspace and solves the corresponding problem.
   *
   * @param parameters new parameters for the R space
   * @param action corresponding action
   * @param pparams P space components of the working set solutions
   * @return
   */
  size_t add_vector(const VecRef<R>& parameters, const VecRef<R>& actions, std::vector<VectorP>& pparams,
                    fapply_on_p_type& apply_p) {
    m_logger->msg("IterativeSolverTemplate::add_vector  iteration = " + std::to_string(m_stats->iterations) +
                      ", apply_p = " + std::to_string(bool(apply_p)),
                  Logger::Trace);
    m_logger->msg("IterativeSolverTemplate::add_vector  size of {params, actions, working_set} = " +
                      std::to_string(parameters.size()) + ", " + std::to_string(actions.size()) + ", " +
                      std::to_string(m_working_set.size()) + ", ",
                  Logger::Debug);
    auto nW = std::min(m_working_set.size(), parameters.size());
    auto cwparams = cwrap<R>(begin(parameters), begin(parameters) + nW);
    auto cwactions = cwrap<R>(begin(actions), begin(actions) + nW);
    m_stats->r_creations += nW;
    m_xspace->update_qspace(cwparams, cwactions);
    return solve_and_generate_working_set(parameters, actions, pparams, apply_p);
  }

public:
  size_t add_vector(const VecRef<R>& parameters, const VecRef<R>& action, fapply_on_p_type& apply_p) override {
    auto pparams = std::vector<VectorP>{};
    return add_vector(parameters, action, pparams, apply_p);
  }

  size_t add_vector(const VecRef<R>& parameters, const VecRef<R>& action, std::vector<VectorP>& pparams) override {
    auto apply_p = fapply_on_p_type{};
    return add_vector(parameters, action, pparams, apply_p);
  }

  size_t add_vector(const VecRef<R>& parameters, const VecRef<R>& action) override {
    auto pparams = std::vector<VectorP>{};
    auto apply_p = fapply_on_p_type{};
    return add_vector(parameters, action, pparams, apply_p);
  }

  size_t add_vector(std::vector<R>& parameters, std::vector<R>& actions, fapply_on_p_type& apply_p) override {
    return add_vector(wrap(parameters), wrap(actions), apply_p);
  }
  size_t add_vector(std::vector<R>& parameters, std::vector<R>& actions, std::vector<VectorP>& pparams) override {
    return add_vector(wrap(parameters), wrap(actions), pparams);
  }
  size_t add_vector(std::vector<R>& parameters, std::vector<R>& actions) override {
    return add_vector(wrap(parameters), wrap(actions));
  }
  size_t add_vector(R& parameters, R& actions) override {
    auto wparams = std::vector<std::reference_wrapper<R>>{std::ref(parameters)};
    auto wactions = std::vector<std::reference_wrapper<R>>{std::ref(parameters)};
    return add_vector(wparams, wactions);
  }

  // FIXME Currently only works if called on an empty subspace. Either enforce it or generalise.
  size_t add_p(const CVecRef<P>& pparams, const array::Span<value_type>& pp_action_matrix, const VecRef<R>& parameters,
               const VecRef<R>& actions, std::vector<VectorP>& parametersP,
               const fapply_on_p_type& apply_p = fapply_on_p_type{}) override {
    m_xspace->update_pspace(pparams, pp_action_matrix);
    return solve_and_generate_working_set(parameters, actions, parametersP, apply_p);
  };

  void clearP() override {}

  void solution(const std::vector<int>& roots, const VecRef<R>& parameters, const VecRef<R>& residual) override {
    detail::construct_solution(parameters, roots, m_subspace_solver->solutions(), m_xspace->paramsp(),
                               m_xspace->paramsq(), m_xspace->paramsd(), m_xspace->dimensions().oP,
                               m_xspace->dimensions().oQ, m_xspace->dimensions().oD, *m_handlers);
    detail::construct_solution(residual, roots, m_subspace_solver->solutions(), m_xspace->actionsp(),
                               m_xspace->actionsq(), m_xspace->actionsd(), m_xspace->dimensions().oP,
                               m_xspace->dimensions().oQ, m_xspace->dimensions().oD, *m_handlers);
    detail::construct_residual(roots, m_subspace_solver->eigenvalues(), cwrap(parameters), residual, m_handlers->rr());
  };

  void solution(const std::vector<int>& roots, std::vector<R>& parameters, std::vector<R>& residual) override {
    return solution(roots, wrap(parameters), wrap(residual));
  }

  void solution_params(const std::vector<int>& roots, std::vector<R>& parameters) override {
    return solution_params(roots, wrap(parameters));
  }

  void solution_params(const std::vector<int>& roots, const VecRef<R>& parameters) override {
    detail::construct_solution(parameters, roots, m_subspace_solver->solutions(), m_xspace->paramsp(),
                               m_xspace->paramsq(), m_xspace->paramsd(), m_xspace->dimensions().oP,
                               m_xspace->dimensions().oQ, m_xspace->dimensions().oD, *m_handlers);
  };

  // TODO Implement this
  std::vector<size_t> suggest_p(const CVecRef<R>& solution, const CVecRef<R>& residual, size_t max_number,
                                double threshold) override {
    return {};
  }

  const std::vector<int>& working_set() const override { return m_working_set; }

  size_t n_roots() const override { return m_nroots; }

  void set_n_roots(size_t roots) override {
    m_nroots = roots;
    m_working_set.resize(roots);
    std::iota(begin(m_working_set), end(m_working_set), (int)0);
  }

  const std::vector<scalar_type>& errors() const override { return m_errors; }

  const Statistics& statistics() const override { return *m_stats; }

  void report(std::ostream& cout) const override {
    cout << "iteration " << m_stats->iterations;
    if (not m_errors.empty()) {
      auto it_max_error = std::max_element(m_errors.cbegin(), m_errors.cend());
      if (n_roots() > 1)
        cout << ", error[" << std::distance(m_errors.cbegin(), it_max_error) << "] = ";
      else
        cout << ", error = ";
      cout << *it_max_error << std::endl;
    }
  }

  void report() const override { report(std::cout); }

  void set_convergence_threshold(double thresh) override { m_convergence_threshold = thresh; }
  double convergence_threshold() const override { return m_convergence_threshold; }
  //! Access dimensions of the subspace
  const subspace::Dimensions& dimensions() const override { return m_xspace->dimensions(); }

protected:
  IterativeSolverTemplate(std::shared_ptr<subspace::XSpaceI<R, Q, P>> xspace,
                          std::shared_ptr<subspace::SubspaceSolverI<R, Q, P>> solver,
                          std::shared_ptr<ArrayHandlers<R, Q, P>> handlers, std::shared_ptr<Statistics> stats,
                          std::shared_ptr<Logger> logger)
      : m_handlers(std::move(handlers)), m_xspace(std::move(xspace)), m_subspace_solver(std::move(solver)),
        m_stats(std::move(stats)), m_logger(std::move(logger)) {}

  size_t solve_and_generate_working_set(const VecRef<R>& parameters, const VecRef<R>& action,
                                        std::vector<VectorP>& pparams, const fapply_on_p_type& apply_p) {
    m_subspace_solver->solve(*m_xspace, n_roots());
    auto nsol = m_subspace_solver->size();
    std::vector<std::pair<Q, Q>> temp_solutions{};
    for (const auto& batch : detail::parameter_batches(nsol, parameters.size())) {
      auto [start_sol, end_sol] = batch;
      auto roots = std::vector<int>(end_sol - start_sol);
      std::iota(begin(roots), end(roots), start_sol);
      detail::construct_solution(parameters, roots, m_subspace_solver->solutions(), m_xspace->paramsp(),
                                 m_xspace->paramsq(), m_xspace->paramsd(), m_xspace->dimensions().oP,
                                 m_xspace->dimensions().oQ, m_xspace->dimensions().oD, *m_handlers);
      detail::construct_solution(action, roots, m_subspace_solver->solutions(), m_xspace->actionsp(),
                                 m_xspace->actionsq(), m_xspace->actionsd(), m_xspace->dimensions().oP,
                                 m_xspace->dimensions().oQ, m_xspace->dimensions().oD, *m_handlers);
      auto pvectors = detail::construct_vectorP(roots, m_subspace_solver->solutions(), m_xspace->dimensions().oP,
                                                m_xspace->dimensions().nP);
      detail::normalise(roots.size(), parameters, action, m_handlers->rr(), *m_logger);
      if (apply_p) {
        apply_p(pvectors, m_xspace->cparamsp(), action);
      } else {
        detail::remove_p_component(parameters, roots, m_subspace_solver->solutions(), m_xspace->cparamsp(),
                                   m_xspace->dimensions().oP, m_handlers->rp());
        detail::remove_p_component(action, roots, m_subspace_solver->solutions(), m_xspace->cactionsp(),
                                   m_xspace->dimensions().oP, m_handlers->rp());
      }
      detail::construct_residual(roots, m_subspace_solver->eigenvalues(), cwrap(parameters), action, m_handlers->rr());
      auto errors = std::vector<scalar_type>(roots.size(), 0);
      detail::update_errors(errors, cwrap(action), m_handlers->rr());
      for (size_t i = 0; i < roots.size(); ++i)
        temp_solutions.emplace_back(m_handlers->qr().copy(parameters[i]), m_handlers->qr().copy(action[i]));
      m_subspace_solver->set_error(roots, errors);
    }
    m_errors = m_subspace_solver->errors();
    m_working_set = detail::select_working_set(parameters.size(), m_errors, m_convergence_threshold);
    for (size_t i = 0; i < m_working_set.size(); ++i) {
      auto root = m_working_set[i];
      m_handlers->rq().copy(parameters[i], temp_solutions.at(root).first);
      m_handlers->rq().copy(action[i], temp_solutions.at(root).second);
    }
    pparams = detail::construct_vectorP(m_working_set, m_subspace_solver->solutions(), m_xspace->dimensions().oP,
                                        m_xspace->dimensions().nP);
    m_logger->msg("add_vector::errors = ", begin(m_errors), end(m_errors), Logger::Trace, 6);
    return m_working_set.size();
  }

  std::shared_ptr<ArrayHandlers<R, Q, P>> m_handlers;
  std::shared_ptr<subspace::XSpaceI<R, Q, P>> m_xspace;
  std::shared_ptr<subspace::SubspaceSolverI<R, Q, P>> m_subspace_solver;
  std::vector<double> m_errors;
  std::vector<int> m_working_set;
  size_t m_nroots{0};
  double m_convergence_threshold{1.0e-10}; //!< errors less than this mark a converged solution
  std::shared_ptr<Statistics> m_stats;
  std::shared_ptr<Logger> m_logger;
};

} // namespace molpro::linalg::itsolv

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_ITERATIVESOLVERTEMPLATE_H