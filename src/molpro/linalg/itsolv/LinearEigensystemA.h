#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_LINEAREIGENSYSTEMA_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_LINEAREIGENSYSTEMA_H
#include <iterator>
#include <molpro/linalg/itsolv/DSpaceResetter.h>
#include <molpro/linalg/itsolv/IterativeSolverTemplate.h>
#include <molpro/linalg/itsolv/Logger.h>
#include <molpro/linalg/itsolv/propose_rspace.h>
#include <molpro/linalg/itsolv/subspace/SubspaceSolverLinEig.h>
#include <molpro/linalg/itsolv/subspace/XSpace.h>

namespace molpro::linalg::itsolv {
// namespace detail {
// template <class R>
// using SubspaceSolver = subspace::SubspaceSolverLinEig<typename array::ArrayHandler<R, R>::value_type,
//                                                      typename array::ArrayHandler<R, R>::value_type_abs>;
//}

/*!
 * @brief One instance of LinearEigensystem (codename A)
 *
 * This is the leaf and it selects which policies to use.
 *
 * It should implement any left-over methods.
 *
 * @tparam R
 * @tparam Q
 * @tparam P
 */
template <class R, class Q, class P>
class LinearEigensystemA : public IterativeSolverTemplate<LinearEigensystem<R, Q, P>, subspace::XSpace<R, Q, P>> {
public:
  using SolverTemplate = IterativeSolverTemplate<LinearEigensystem<R, Q, P>, subspace::XSpace<R, Q, P>>;
  using typename SolverTemplate::scalar_type;

  explicit LinearEigensystemA(const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers,
                              const std::shared_ptr<Logger>& logger_ = std::make_shared<Logger>())
      : SolverTemplate(subspace::XSpace<R, Q, P>(handlers, logger_),
                       std::static_pointer_cast<subspace::SubspaceSolverI<R, Q, P>>(
                           std::make_shared<subspace::SubspaceSolverLinEig<R, Q, P>>(logger_)),
                       handlers, std::make_shared<Statistics>(), logger_),
        logger(logger_) {
    auto* subspace_solver_lin_eig =
        static_cast<subspace::SubspaceSolverLinEig<R, Q, P>*>(this->m_subspace_solver.get());
    subspace_solver_lin_eig->m_norm_stability_threshold = propose_rspace_norm_thresh;
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
  size_t end_iteration(std::vector<R>& parameters, std::vector<R>& action) override {
    if (m_dspace_resetter.do_reset(this->m_stats->iterations, this->m_xspace.dimensions())) {
      this->m_working_set = m_dspace_resetter.run(parameters, this->m_xspace, this->m_subspace_solver->solutions(),
                                                  *this->m_handlers, *this->m_logger);
    } else {
      this->m_working_set =
          detail::propose_rspace(*static_cast<LinearEigensystem<R, Q, P>*>(this), parameters, action, this->m_xspace,
                                 this->m_subspace_solver->solutions(), *this->m_handlers, *this->m_logger,
                                 propose_rspace_norm_thresh, max_size_qspace);
    }
    this->m_stats->iterations++;
    return this->working_set().size();
  }

  //! Applies the Davidson preconditioner
  void precondition(std::vector<R>& parameters, std::vector<R>& action) const {}

  void set_convergence_threshold(double threshold) { this->m_convergence_threshold = threshold; }

  std::vector<scalar_type> eigenvalues() const override { return this->m_subspace_solver->eigenvalues(); }

  std::vector<scalar_type> working_set_eigenvalues() const {
    auto eval = std::vector<scalar_type>{};
    for (auto i : this->working_set()) {
      eval.emplace_back(this->m_subspace_solver->eigenvalues().at(i));
    }
    return eval;
  }

  void report() const override {
    SolverTemplate::report();
    molpro::cout << "errors " << std::scientific;
    auto& err = this->m_errors;
    std::copy(begin(err), end(err), std::ostream_iterator<scalar_type>(molpro::cout, ", "));
    molpro::cout << std::endl;
    molpro::cout << "eigenvalues ";
    auto ev = eigenvalues();
    molpro::cout << std::fixed << std::setprecision(14);
    std::copy(begin(ev), end(ev), std::ostream_iterator<scalar_type>(molpro::cout, ", "));
    molpro::cout << std::defaultfloat << std::endl;
  }

  //! Set the period in iterations for resetting the D space
  void set_reset_D(size_t n) { m_dspace_resetter.set_nreset(n); }
  //! Set the maximum size of Q space after resetting the D space
  void set_reset_D_maxQ_size(size_t n) { m_dspace_resetter.set_max_Qsize(n); }

  std::shared_ptr<Logger> logger;
  double propose_rspace_norm_thresh = 1e-6; //!< vectors with norm less than threshold can be considered null
  unsigned int max_size_qspace = std::numeric_limits<unsigned int>::max(); //!< maximum size of Q space
protected:
  detail::DSpaceResetter<Q> m_dspace_resetter; //!< resets D space
};

} // namespace molpro::linalg::itsolv

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_LINEAREIGENSYSTEMA_H