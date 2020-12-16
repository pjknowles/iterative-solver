#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_ITERATIVESOLVERTEMPLATE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_ITERATIVESOLVERTEMPLATE_H
#include <iostream>
#include <molpro/iostream.h>
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/Logger.h>
#include <molpro/linalg/itsolv/subspace/ISubspaceSolver.h>
#include <molpro/linalg/itsolv/subspace/IXSpace.h>
#include <molpro/linalg/itsolv/subspace/Matrix.h>
#include <molpro/linalg/itsolv/subspace/util.h>
#include <molpro/linalg/itsolv/wrap.h>
#include <stack>

namespace molpro::linalg::itsolv {
namespace detail {

inline std::vector<std::pair<size_t, size_t>> parameter_batches(const size_t nsol, const size_t nparam) {
  auto batches = std::vector<std::pair<size_t, size_t>>{};
  if (nparam && nsol) {
    auto n_batch = nsol / nparam + (nsol % nparam ? 1 : 0);
    for (size_t ib = 0, start_sol = 0, end_sol = 0; ib < n_batch; ++ib, start_sol = end_sol) {
      end_sol = std::min(start_sol + nparam, nsol);
      batches.emplace_back(start_sol, end_sol);
    }
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
void update_errors(std::vector<T>& errors, const CVecRef<R>& residual, array::ArrayHandler<R, R>& handler) {
  assert(residual.size() >= errors.size());
  for (size_t i = 0; i < errors.size(); ++i) {
    auto a = handler.dot(residual[i], residual[i]);
    errors[i] = std::sqrt(std::abs(a));
  }
}

template <typename T>
std::vector<int> select_working_set(const size_t nw, const std::vector<T>& errors, const T threshold,
                                    const std::vector<T>& value_errors, const T value_threshold) {
  auto ordered_errors = std::multimap<T, size_t, std::greater<T>>{};
  for (size_t i = 0; i < errors.size(); ++i) {
    if (errors[i] > threshold or (i < value_errors.size() and value_errors[i] > value_threshold))
      ordered_errors.emplace(errors[i], i);
  }
  auto working_set = std::vector<int>{};
  auto end = (ordered_errors.size() < nw ? ordered_errors.end() : next(begin(ordered_errors), nw));
  std::transform(begin(ordered_errors), end, std::back_inserter(working_set), [](const auto& el) { return el.second; });
  return working_set;
}

} // namespace detail

/*!
 * @brief Implements IterativeSolver interface that is common to all solvers
 *
 * @tparam Solver one of the solvers that inherits from IterativeSolver.
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

public:
  size_t add_vector(const VecRef<R>& parameters, const VecRef<R>& actions) override {
    m_logger->msg("IterativeSolverTemplate::add_vector  iteration = " + std::to_string(m_stats->iterations),
                  Logger::Trace);
    m_logger->msg("IterativeSolverTemplate::add_vector  size of {params, actions, working_set} = " +
                      std::to_string(parameters.size()) + ", " + std::to_string(actions.size()) + ", " +
                      std::to_string(m_working_set.size()) + ", ",
                  Logger::Debug);
    if (m_xspace->dimensions().nP != 0 && !m_apply_p)
      throw std::runtime_error(
          "Solver contains P space but no valid apply_p function. Make sure add_p was called correctly.");
    auto nW = std::min(m_working_set.size(), parameters.size());
    auto cwparams = cwrap<R>(begin(parameters), begin(parameters) + nW);
    auto cwactions = cwrap<R>(begin(actions), begin(actions) + nW);
    m_stats->r_creations += nW;
    m_xspace->update_qspace(cwparams, cwactions);
    return solve_and_generate_working_set(parameters, actions);
  }

public:
  size_t add_vector(std::vector<R>& parameters, std::vector<R>& actions) override {
    return add_vector(wrap(parameters), wrap(actions));
  }
  size_t add_vector(R& parameters, R& actions) override {
    auto wparams = std::vector<std::reference_wrapper<R>>{std::ref(parameters)};
    auto wactions = std::vector<std::reference_wrapper<R>>{std::ref(actions)};
    return add_vector(wparams, wactions);
  }

  // FIXME Currently only works if called on an empty subspace. Either enforce it or generalise.
  size_t add_p(const CVecRef<P>& pparams, const array::Span<value_type>& pp_action_matrix, const VecRef<R>& parameters,
               const VecRef<R>& actions, fapply_on_p_type apply_p) override {
    if (not pparams.empty() and pparams.size() < n_roots())
      throw std::runtime_error("P space must be empty or at least as large as number of roots sought");
    if (apply_p)
      m_apply_p = std::move(apply_p);
    m_xspace->update_pspace(pparams, pp_action_matrix);
    return solve_and_generate_working_set(parameters, actions);
  };

  void clearP() override {}

  void solution(const std::vector<int>& roots, const VecRef<R>& parameters, const VecRef<R>& residual) override {
    check_consistent_number_of_roots_and_solutions(roots, parameters.size());
    detail::construct_solution(parameters, roots, m_subspace_solver->solutions(), m_xspace->paramsp(),
                               m_xspace->paramsq(), m_xspace->paramsd(), m_xspace->dimensions().oP,
                               m_xspace->dimensions().oQ, m_xspace->dimensions().oD, *m_handlers);
    detail::construct_solution(residual, roots, m_subspace_solver->solutions(), {}, m_xspace->actionsq(),
                               m_xspace->actionsd(), m_xspace->dimensions().oP, m_xspace->dimensions().oQ,
                               m_xspace->dimensions().oD, *m_handlers);
    auto pvectors = detail::construct_vectorP(roots, m_subspace_solver->solutions(), m_xspace->dimensions().oP,
                                              m_xspace->dimensions().nP);
    if (m_normalise_solution)
      detail::normalise(roots.size(), parameters, residual, m_handlers->rr(), *m_logger);
    if (m_apply_p)
      m_apply_p(pvectors, m_xspace->cparamsp(), residual);
    construct_residual(roots, cwrap(parameters), residual);
  };

  void solution(const std::vector<int>& roots, std::vector<R>& parameters, std::vector<R>& residual) override {
    return solution(roots, wrap(parameters), wrap(residual));
  }

  void solution_params(const std::vector<int>& roots, std::vector<R>& parameters) override {
    return solution_params(roots, wrap(parameters));
  }

  void solution_params(const std::vector<int>& roots, const VecRef<R>& parameters) override {
    check_consistent_number_of_roots_and_solutions(roots, parameters.size());
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

  void set_options(const Options& options) override {
    if (options.n_roots)
      set_n_roots(options.n_roots.value());
    if (options.convergence_threshold)
      set_convergence_threshold(options.convergence_threshold.value());
  }

  std::shared_ptr<Options> get_options() const override {
    auto options = std::make_shared<Options>();
    options->n_roots = n_roots();
    options->convergence_threshold = convergence_threshold();
    return options;
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

  void report() const override { report(molpro::cout); }

  void set_convergence_threshold(double thresh) override { m_convergence_threshold = thresh; }
  double convergence_threshold() const override { return m_convergence_threshold; }
  void set_convergence_threshold_value(double thresh) override { m_convergence_threshold_value = thresh; }
  double convergence_threshold_value() const override { return m_convergence_threshold_value; }
  //! Access dimensions of the subspace
  const subspace::Dimensions& dimensions() const override { return m_xspace->dimensions(); }

protected:
  IterativeSolverTemplate(std::shared_ptr<subspace::IXSpace<R, Q, P>> xspace,
                          std::shared_ptr<subspace::ISubspaceSolver<R, Q, P>> solver,
                          std::shared_ptr<ArrayHandlers<R, Q, P>> handlers, std::shared_ptr<Statistics> stats,
                          std::shared_ptr<Logger> logger)
      : m_handlers(std::move(handlers)), m_xspace(std::move(xspace)), m_subspace_solver(std::move(solver)),
        m_stats(std::move(stats)), m_logger(std::move(logger)) {}

  //! Implementation class should overload this to set errors in the current values (e.g. change in eigenvalues)
  virtual void set_value_errors() {}
  //! Constructs residual for given roots provided their parameters and actions
  virtual void construct_residual(const std::vector<int>& roots, const CVecRef<R>& params,
                                  const VecRef<R>& actions) = 0;

  /*!
   * @brief Solves the subspace problems and selects the working set of roots, returning their parameters and residual
   * in parameters and action
   * @param parameters container for storing parameters of the working set
   * @param action container for storing the residual of the working set
   * @param apply_p function that accumulates action from the P space projection of parameters
   * @return size of the working set
   */
  size_t solve_and_generate_working_set(const VecRef<R>& parameters, const VecRef<R>& action) {
    m_subspace_solver->solve(*m_xspace, n_roots());
    auto nsol = m_subspace_solver->size();
    std::vector<std::pair<Q, Q>> temp_solutions{};
    for (const auto& batch : detail::parameter_batches(nsol, parameters.size())) {
      auto [start_sol, end_sol] = batch;
      auto roots = std::vector<int>(end_sol - start_sol);
      std::iota(begin(roots), end(roots), start_sol);
      solution(roots, parameters, action);
      auto errors = std::vector<scalar_type>(roots.size(), 0);
      detail::update_errors(errors, cwrap(action), m_handlers->rr());
      for (size_t i = 0; i < roots.size(); ++i)
        temp_solutions.emplace_back(m_handlers->qr().copy(parameters[i]), m_handlers->qr().copy(action[i]));
      m_subspace_solver->set_error(roots, errors);
    }
    set_value_errors();
    m_errors = m_subspace_solver->errors();
    m_working_set = detail::select_working_set(parameters.size(), m_errors, m_convergence_threshold, m_value_errors,
                                               m_convergence_threshold_value);
    for (size_t i = 0; i < m_working_set.size(); ++i) {
      auto root = m_working_set[i];
      m_handlers->rq().copy(parameters[i], temp_solutions.at(root).first);
      m_handlers->rq().copy(action[i], temp_solutions.at(root).second);
    }
    m_logger->msg("add_vector::errors = ", begin(m_errors), end(m_errors), Logger::Trace, 6);
    return m_working_set.size();
  }

  template <typename I>
  void check_consistent_number_of_roots_and_solutions(const std::vector<I>& roots, const size_t nparams) {
    if (roots.size() > nparams)
      throw std::runtime_error("asking for more roots than parameters");
    if (!roots.empty() && *std::max_element(roots.begin(), roots.end()) >= m_subspace_solver->solutions().size())
      throw std::runtime_error("asking for more roots than there are solutions");
  }

  std::shared_ptr<ArrayHandlers<R, Q, P>> m_handlers;   //!< Array handlers
  std::shared_ptr<subspace::IXSpace<R, Q, P>> m_xspace; //!< manages the subspace and associated data
  std::stack<scalar_type> m_values;                     //! function values in this and previous iterations
  std::shared_ptr<subspace::ISubspaceSolver<R, Q, P>> m_subspace_solver; //!< solves the subspace problem
  std::vector<double> m_errors;                                          //!< errors from the most recent solution
  std::vector<double> m_value_errors;                                    //!< value errors from the most recent solution
  std::vector<int> m_working_set;                                        //!< indices of roots in the working set
  size_t m_nroots{0};                      //!< number of roots the solver is searching for
  double m_convergence_threshold{1.0e-10}; //!< residual norms less than this mark a converged solution
  double m_convergence_threshold_value{
      std::numeric_limits<double>::max()}; //!< value changes less than this mark a converged solution
  std::shared_ptr<Statistics> m_stats;     //!< accumulates statistics of operations performed by the solver
  std::shared_ptr<Logger> m_logger;        //!< logger
  bool m_normalise_solution = false;       //!< whether to normalise the solutions
  fapply_on_p_type m_apply_p = {};         //!< function that evaluates effect of action on the P space projection
};

} // namespace molpro::linalg::itsolv

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_ITERATIVESOLVERTEMPLATE_H
