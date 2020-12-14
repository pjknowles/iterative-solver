#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_LINEAREIGENSYSTEM_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_LINEAREIGENSYSTEM_H
#include <iterator>
#include <molpro/linalg/itsolv/CastOptions.h>
#include <molpro/linalg/itsolv/DSpaceResetter.h>
#include <molpro/linalg/itsolv/IterativeSolverTemplate.h>
#include <molpro/linalg/itsolv/Logger.h>
#include <molpro/linalg/itsolv/propose_rspace.h>
#include <molpro/linalg/itsolv/subspace/SubspaceSolverLinEig.h>
#include <molpro/linalg/itsolv/subspace/XSpace.h>

namespace molpro::linalg::itsolv {

/*!
 * @brief One specific implementation of LinearEigensystem (codename A)
 *
 * This implementation uses our own D-space algorithm.
 *
 * TODO add more documentation and examples
 *
 * @tparam R
 * @tparam Q
 * @tparam P
 */
template <class R, class Q, class P>
class LinearEigensystem : public IterativeSolverTemplate<ILinearEigensystem, R, Q, P> {
public:
  using SolverTemplate = IterativeSolverTemplate<ILinearEigensystem, R, Q, P>;
  using typename SolverTemplate::scalar_type;
  using IterativeSolverTemplate<ILinearEigensystem, R, Q, P>::report;

  explicit LinearEigensystem(const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers,
                             const std::shared_ptr<Logger>& logger_ = std::make_shared<Logger>())
      : SolverTemplate(std::make_shared<subspace::XSpace<R, Q, P>>(handlers, logger_),
                       std::static_pointer_cast<subspace::ISubspaceSolver<R, Q, P>>(
                           std::make_shared<subspace::SubspaceSolverLinEig<R, Q, P>>(logger_)),
                       handlers, std::make_shared<Statistics>(), logger_),
        logger(logger_) {
    set_hermiticity(m_hermiticity);
    this->m_normalise_solution = false;
  }

  /*!
   * \brief Proposes new parameters for the subspace from the preconditioned residuals.
   *
   * After add_vector solves the subspace problem it returns the new solution and residual. The residual should be
   * preconditioned, i.e. using Davidson method, to accelerate convergence. The updated residual is used in this
   * function to propose new Q space parameters orthonormal to the old space. They are returned in parameters so that
   * corresponding actions can be calculated and used in add_vector in the next iteration.
   *
   * Every n_reset_D iterations the D space has to be reset. If working set does not cover all solutions than resetting
   * is done over multiple iterations. During this time, there are no deletions from the Q space.
   * This is done by temporarily increasing maximum allowed size of Q space.
   *
   * @param parameters output new parameters for the subspace.
   * @param residual preconditioned residuals.
   * @return number of significant parameters to calculate the action for
   */
  size_t end_iteration(const VecRef<R>& parameters, const VecRef<R>& action) override {
    if (m_dspace_resetter.do_reset(this->m_stats->iterations, this->m_xspace->dimensions())) {
      m_resetting_in_progress = true;
      this->m_working_set = m_dspace_resetter.run(parameters, *this->m_xspace, this->m_subspace_solver->solutions(),
                                                  propose_rspace_norm_thresh, propose_rspace_svd_thresh,
                                                  *this->m_handlers, *this->m_logger);
    } else {
      m_resetting_in_progress = false;
      this->m_working_set = detail::propose_rspace(*this, parameters, action, *this->m_xspace, *this->m_subspace_solver,
                                                   *this->m_handlers, *this->m_logger, propose_rspace_svd_thresh,
                                                   propose_rspace_norm_thresh, m_max_size_qspace);
    }
    this->m_stats->iterations++;
    return this->working_set().size();
  }
  size_t end_iteration(std::vector<R>& parameters, std::vector<R>& action) override {
    return end_iteration(wrap(parameters), wrap(action));
  }

  //! Applies the Davidson preconditioner
  void precondition(std::vector<R>& parameters, std::vector<R>& action) const {}

  std::vector<scalar_type> eigenvalues() const override { return this->m_subspace_solver->eigenvalues(); }

  std::vector<scalar_type> working_set_eigenvalues() const override {
    auto eval = std::vector<scalar_type>{};
    for (auto i : this->working_set()) {
      eval.emplace_back(this->m_subspace_solver->eigenvalues().at(i));
    }
    return eval;
  }

  void set_value_errors() override {
    auto current_values = this->m_subspace_solver->eigenvalues();
    this->m_value_errors.assign(current_values.size(), std::numeric_limits<scalar_type>::max());
    for (size_t i = 0; i < std::min(m_last_values.size(), current_values.size()); i++)
      this->m_value_errors[i] = std::abs(current_values[i] - m_last_values[i]);
    if (!m_resetting_in_progress)
      m_last_values = current_values;
  }

  void report(std::ostream& cout) const override {
    SolverTemplate::report(cout);
    cout << "errors " << std::scientific;
    auto& err = this->m_errors;
    std::copy(begin(err), end(err), std::ostream_iterator<scalar_type>(molpro::cout, ", "));
    cout << std::endl;
    cout << "eigenvalues ";
    auto ev = eigenvalues();
    cout << std::fixed << std::setprecision(14);
    std::copy(begin(ev), end(ev), std::ostream_iterator<scalar_type>(molpro::cout, ", "));
    cout << std::defaultfloat << std::endl;
  }

  //! Set the period in iterations for resetting the D space
  void set_reset_D(size_t n) { m_dspace_resetter.set_nreset(n); }
  size_t get_reset_D() const { return m_dspace_resetter.get_nreset(); }
  //! Set the maximum size of Q space after resetting the D space
  void set_reset_D_maxQ_size(size_t n) { m_dspace_resetter.set_max_Qsize(n); }
  int get_reset_D_maxQ_size() const { return m_dspace_resetter.get_max_Qsize(); }
  int get_max_size_qspace() const { return m_max_size_qspace; }
  void set_max_size_qspace(int n) {
    m_max_size_qspace = n;
    if (m_dspace_resetter.get_max_Qsize() > m_max_size_qspace)
      m_dspace_resetter.set_max_Qsize(m_max_size_qspace);
  }
  void set_hermiticity(bool hermitian) override {
    m_hermiticity = hermitian;
    auto xspace = std::dynamic_pointer_cast<subspace::XSpace<R, Q, P>>(this->m_xspace);
    xspace->set_hermiticity(hermitian);
    auto subspace_solver = std::dynamic_pointer_cast<subspace::SubspaceSolverLinEig<R, Q, P>>(this->m_subspace_solver);
    subspace_solver->set_hermiticity(hermitian);
  }
  bool get_hermiticity() const override { return m_hermiticity; }

  void set_options(const std::shared_ptr<Options>& options) override {
    SolverTemplate::set_options(options);
    auto opt = CastOptions::LinearEigensystem(options);
    if (opt) {
      if (opt->reset_D)
        set_reset_D(opt->reset_D.value());
      if (opt->reset_D_max_Q_size)
        set_reset_D_maxQ_size(opt->reset_D_max_Q_size.value());
      if (opt->max_size_qspace)
        set_max_size_qspace(opt->max_size_qspace.value());
      if (opt->norm_thresh)
        propose_rspace_norm_thresh = opt->norm_thresh.value();
      if (opt->svd_thresh)
        propose_rspace_svd_thresh = opt->svd_thresh.value();
      if (opt->hermiticity)
        set_hermiticity(opt->hermiticity.value());
    }
  }

  std::shared_ptr<Options> get_options() const override {
    auto opt = std::make_shared<LinearEigensystemOptions>();
    opt->copy(*SolverTemplate::get_options());
    opt->reset_D = get_reset_D();
    opt->reset_D_max_Q_size = get_reset_D_maxQ_size();
    opt->max_size_qspace = get_max_size_qspace();
    opt->norm_thresh = propose_rspace_norm_thresh;
    opt->svd_thresh = propose_rspace_svd_thresh;
    opt->hermiticity = get_hermiticity();
    return opt;
  }

  std::shared_ptr<Logger> logger;
  double propose_rspace_norm_thresh = 1e-10; //!< vectors with norm less than threshold can be considered null.
  double propose_rspace_svd_thresh = 1e-12;  //!< the smallest singular value in the subspace that can be allowed when
                                             //!< constructing the working set. Smaller singular values will lead to
                                             //!< deletion of parameters from the Q space
protected:
  void construct_residual(const std::vector<int>& roots, const CVecRef<R>& params, const VecRef<R>& actions) override {
    assert(params.size() >= roots.size());
    const auto& eigvals = eigenvalues();
    for (size_t i = 0; i < roots.size(); ++i)
      this->m_handlers->rr().axpy(-eigvals.at(roots[i]), params.at(i), actions.at(i));
  }

  int m_max_size_qspace = std::numeric_limits<int>::max(); //!< maximum size of Q space
  detail::DSpaceResetter<Q> m_dspace_resetter;             //!< resets D space
  bool m_hermiticity = false;                              //!< whether the problem is hermitian or not
  std::vector<double> m_last_values;                       //!< The values from the previous iteration
  bool m_resetting_in_progress = false;                    //!< whether D space resetting is in progress
};

} // namespace molpro::linalg::itsolv

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_LINEAREIGENSYSTEM_H