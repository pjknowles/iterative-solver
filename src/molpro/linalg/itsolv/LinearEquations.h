#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_LINEAREQUATIONS_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_LINEAREQUATIONS_H
#include <molpro/linalg/itsolv/DSpaceResetter.h>
#include <molpro/linalg/itsolv/IterativeSolverTemplate.h>
#include <molpro/linalg/itsolv/propose_rspace.h>

namespace molpro::linalg::itsolv {
/*!
 * @brief Solves a system of linear equation, A x = b
 * @tparam R
 * @tparam Q
 * @tparam P
 */
template <class R, class Q, class P>
class LinearEquations : public IterativeSolverTemplate<ILinearEquations, R, Q, P> {
public:
  using SolverTemplate = IterativeSolverTemplate<ILinearEquations, R, Q, P>;
  using SolverTemplate ::report;
  using typename SolverTemplate::scalar_type;

  explicit LinearEquations(const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers,
                           const std::shared_ptr<Logger>& logger_ = std::make_shared<Logger>())
      : SolverTemplate(std::make_shared<subspace::XSpace<R, Q, P>>(handlers, logger_),
                       std::static_pointer_cast<subspace::ISubspaceSolver<R, Q, P>>(
                           std::make_shared<subspace::SubspaceSolverLinEig<R, Q, P>>(logger_)),
                       handlers, std::make_shared<Statistics>(), logger_),
        logger(logger_) {
    set_hermiticity(m_hermiticity);
  }

  size_t end_iteration(const VecRef<R>& parameters, const VecRef<R>& action) override {
    if (m_dspace_resetter.do_reset(this->m_stats->iterations, this->m_xspace->dimensions())) {
      this->m_working_set = m_dspace_resetter.run(parameters, *this->m_xspace, this->m_subspace_solver->solutions(),
                                                  m_norm_thresh, m_svd_thresh, *this->m_handlers, *this->m_logger);
    } else {
      this->m_working_set =
          detail::propose_rspace(*this, parameters, action, *this->m_xspace, *this->m_subspace_solver,
                                 *this->m_handlers, *this->m_logger, m_svd_thresh, m_norm_thresh, m_max_size_qspace);
    }
    this->m_stats->iterations++;
    return this->working_set().size();
  }

  // FIXME move this to the template
  size_t end_iteration(std::vector<R>& parameters, std::vector<R>& action) override {
    return end_iteration(wrap(parameters), wrap(action));
  }

  void add_equations(const CVecRef<R>& rhs) override {
    this->set_n_roots(rhs.size());
    for (const auto& r : rhs) {
      m_rhs.emplace_back(this->m_handlers->qr().copy(r));
    }
  }

  void add_equations(const R& rhs) override { add_equations(cwrap_arg(rhs)); }

  const std::vector<Q>& rhs() const override { return m_rhs; }

  //! Set threshold on the norm of parameters that should be considered null
  void set_norm_thresh(double thresh) { m_norm_thresh = thresh; }
  double get_norm_thresh() const { return m_norm_thresh; }
  //! Set the smallest singular value in the subspace that can be allowed when
  //! constructing the working set. Smaller singular values will lead to deletion of parameters
  void set_svd_thresh(double thresh) { m_svd_thresh = thresh; }
  double get_svd_thresh() const { return m_svd_thresh; }
  //! Set the period in iterations for resetting the D space
  void set_reset_D(size_t n) { m_dspace_resetter.set_nreset(n); }
  size_t get_reset_D() const { return m_dspace_resetter.get_nreset(); }
  //! Set the maximum size of Q space after resetting the D space
  void set_reset_D_maxQ_size(size_t n) { m_dspace_resetter.set_max_Qsize(n); }
  int get_reset_D_maxQ_size() const { return m_dspace_resetter.get_max_Qsize(); }
  //! Set a limit on the maximum size of Q space. This does not include the size of the working space (R) and the D
  //! space
  void set_max_size_qspace(int n) {
    m_max_size_qspace = n;
    if (m_dspace_resetter.get_max_Qsize() > m_max_size_qspace)
      m_dspace_resetter.set_max_Qsize(m_max_size_qspace);
  }
  int get_max_size_qspace() const { return m_max_size_qspace; }
  void set_hermiticity(bool hermitian) override {
    m_hermiticity = hermitian;
    auto xspace = std::dynamic_pointer_cast<subspace::XSpace<R, Q, P>>(this->m_xspace);
    xspace->set_hermiticity(hermitian);
    auto subspace_solver = std::dynamic_pointer_cast<subspace::SubspaceSolverLinEig<R, Q, P>>(this->m_subspace_solver);
    subspace_solver->set_hermiticity(hermitian);
  }
  bool get_hermiticity() const override { return m_hermiticity; }

  std::shared_ptr<Logger> logger;

protected:
  void construct_residual(const std::vector<int>& roots, const CVecRef<R>& params, const VecRef<R>& actions) override {
    assert(params.size() >= roots.size());
    for (size_t i = 0; i < roots.size(); ++i)
      this->m_handlers->rq().axpy(-1, m_rhs.at(roots[i]), actions.at(i));
  }

  std::vector<Q> m_rhs;         //!< Right hand side vectors
  double m_norm_thresh = 1e-10; //!< vectors with norm less than threshold can be considered null.
  double m_svd_thresh = 1e-8;   //!< svd values smaller than this mark the null space
  int m_max_size_qspace = std::numeric_limits<int>::max(); //!< maximum size of Q space
  detail::DSpaceResetter<Q> m_dspace_resetter;             //!< resets D space
  bool m_hermiticity = true;                               //!< whether the problem is hermitian or not
};

} // namespace molpro::linalg::itsolv
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_LINEAREQUATIONS_H
