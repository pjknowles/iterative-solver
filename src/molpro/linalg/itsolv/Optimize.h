#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_OPTIMIZE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_OPTIMIZE_H
#include <molpro/linalg/itsolv/CastOptions.h>
#include <molpro/linalg/itsolv/DSpaceResetter.h>
#include <molpro/linalg/itsolv/IterativeSolverTemplate.h>
#include <molpro/linalg/itsolv/propose_rspace.h>
#include <molpro/linalg/itsolv/subspace/XSpace.h>
#include <molpro/linalg/itsolv/subspace/SubspaceSolverOpt.h>

namespace molpro::linalg::itsolv {
/*!
 * @brief A class that optimises a function using a Quasi-Newton or other method.
 *
 * @todo Explain some of the theory
 *
 * @tparam R The class encapsulating solution and residual vectors
 * @tparam Q Used internally as a class for storing vectors on backing store
 */
template <class R, class Q>
class Optimize : public IterativeSolverTemplate<IOptimize, R, Q, R> {
public:
  using P = R;
  using SolverTemplate = IterativeSolverTemplate<IOptimize, R, Q, P>;
  using SolverTemplate ::report;
  using typename SolverTemplate::value_type;
  using typename SolverTemplate::value_type_abs;

  explicit Optimize(const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers,
                           const std::shared_ptr<Logger>& logger_ = std::make_shared<Logger>())
      : SolverTemplate(std::make_shared<subspace::XSpace<R, Q, P>>(handlers, logger_),
                       std::static_pointer_cast<subspace::ISubspaceSolver<R, Q, P>>(
                           std::make_shared<subspace::SubspaceSolverOpt<R, Q, P>>(logger_)),
                       handlers, std::make_shared<Statistics>(), logger_),
        logger(logger_) {
    set_hermiticity(m_hermiticity);
    this->m_normalise_solution = false;
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

  //! Set threshold on the norm of parameters that should be considered null
  void set_norm_thresh(double thresh) { m_norm_thresh = thresh; }
  double get_norm_thresh() const { return m_norm_thresh; }
  //! Set the smallest singular value in the subspace that can be allowed when
  //! constructing the working set. Smaller singular values will lead to deletion of parameters
  void set_svd_thresh(double thresh) { m_svd_thresh = thresh; }
  double get_svd_thresh() const { return m_svd_thresh; }
  //! Set a limit on the maximum size of Q space. This does not include the size of the working space (R) and the D
  //! space
  void set_max_size_qspace(int n) {
    m_max_size_qspace = n;
  }
  int get_max_size_qspace() const { return m_max_size_qspace; }

  void set_options(const std::shared_ptr<Options>& options) override {
    SolverTemplate::set_options(options);
    auto opt = CastOptions::Optimize(options);
    if (opt) {
      if (opt->reset_D)
        set_reset_D(opt->reset_D.value());
      if (opt->reset_D_max_Q_size)
        set_reset_D_maxQ_size(opt->reset_D_max_Q_size.value());
      if (opt->max_size_qspace)
        set_max_size_qspace(opt->max_size_qspace.value());
      if (opt->norm_thresh)
        set_norm_thresh(opt->norm_thresh.value());
      if (opt->svd_thresh)
        set_svd_thresh(opt->svd_thresh.value());
      if (opt->hermiticity)
        set_hermiticity(opt->hermiticity.value());
      if (opt->augmented_hessian)
        set_augmented_hessian(opt->augmented_hessian.value());
    }
  }

  std::shared_ptr<Options> get_options() const override {
    auto opt = std::make_shared<OptimizeOptions>();
    opt->copy(*SolverTemplate::get_options());
    opt->reset_D = get_reset_D();
    opt->reset_D_max_Q_size = get_reset_D_maxQ_size();
    opt->max_size_qspace = get_max_size_qspace();
    opt->norm_thresh = get_norm_thresh();
    opt->svd_thresh = get_svd_thresh();
    opt->hermiticity = get_hermiticity();
    opt->augmented_hessian = get_augmented_hessian();
    return opt;
  }

  void report(std::ostream& cout) const override {
    SolverTemplate::report(cout);
    cout << "errors " << std::scientific;
    auto& err = this->m_errors;
    std::copy(begin(err), end(err), std::ostream_iterator<value_type_abs>(molpro::cout, ", "));
    cout << std::defaultfloat << std::endl;
  }
  std::shared_ptr<Logger> logger;

protected:
  // FIXME The scale is fixed by the norm of RHS, but if RHS=0 there is no reference. We could use the norm of params
  void construct_residual(const std::vector<int>& roots, const CVecRef<R>& params, const VecRef<R>& actions) override {
    assert(params.size() >= roots.size());
    const auto& norm = std::dynamic_pointer_cast<subspace::XSpace<R, Q, P>>(this->m_xspace)->rhs_norm();
    for (size_t i = 0; i < roots.size(); ++i) {
      const auto ii = roots[i];
      this->m_handlers->rq().axpy(-1, rhs().at(ii), actions.at(i));
      if (norm.at(ii) != 0) {
        auto scal = 1 / norm[ii];
        this->m_handlers->rr().scal(scal, actions.at(i));
      }
    }
  }

  double m_norm_thresh = 1e-10; //!< vectors with norm less than threshold can be considered null.
  double m_svd_thresh = 1e-12;  //!< svd values smaller than this mark the null space
  int m_max_size_qspace = std::numeric_limits<int>::max(); //!< maximum size of Q space
  detail::DSpaceResetter<Q> m_dspace_resetter;             //!< resets D space
  bool m_hermiticity = true;                               //!< whether the problem is hermitian or not
};

} // namespace molpro::linalg::itsolv
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_OPTIMIZE_H
