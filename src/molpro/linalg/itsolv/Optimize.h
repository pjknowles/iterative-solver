#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_OPTIMIZE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_OPTIMIZE_H
#include <molpro/linalg/itsolv/CastOptions.h>
#include <molpro/linalg/itsolv/DSpaceResetter.h>
#include <molpro/linalg/itsolv/IterativeSolverTemplate.h>
#include <molpro/linalg/itsolv/propose_rspace.h>
#include <molpro/linalg/itsolv/subspace/SubspaceSolverOpt.h>
#include <molpro/linalg/itsolv/subspace/XSpace.h>

namespace molpro::linalg::itsolv {
/*!
 * @brief A class that optimises a function using a Quasi-Newton or other method.
 *
 * @todo Explain some of the theory
 *
 * @tparam R The class encapsulating solution and residual vectors
 * @tparam Q Used internally as a class for storing vectors on backing store
 */
template <class R, class Q, class P = R>
class Optimize : public IterativeSolverTemplate<IOptimize, R, Q, P> {
public:
  using SolverTemplate = IterativeSolverTemplate<IOptimize, R, Q, P>;
  using SolverTemplate ::report;
  using typename SolverTemplate::scalar_type;
  using typename SolverTemplate::value_type;
  using typename SolverTemplate::value_type_abs;

  explicit Optimize(const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers,
                    const std::shared_ptr<Logger>& logger_ = std::make_shared<Logger>())
      : SolverTemplate(std::make_shared<subspace::XSpace<R, Q, P>>(handlers, logger_),
                       std::static_pointer_cast<subspace::ISubspaceSolver<R, Q, P>>(
                           std::make_shared<subspace::SubspaceSolverOpt<R, Q, P>>(logger_)),
                       handlers, std::make_shared<Statistics>(), logger_),
        logger(logger_) {}

  size_t end_iteration(const VecRef<R>& parameters, const VecRef<R>& action) override {
    // TODO implement for other optimiser variants
//    this->m_working_set =
//        detail::propose_rspace(*this, parameters, action, *this->m_xspace, *this->m_subspace_solver, *this->m_handlers,
//                               *this->m_logger, m_svd_thresh, m_norm_thresh, m_max_size_qspace);
    this->solution_params(this->m_working_set,parameters);
    this->m_handlers->rr().axpy(1,action.front(),parameters.front());
    if (this->m_errors.front() < this->m_convergence_threshold)
      this->m_working_set.clear();
    else
      this->m_working_set.assign(1, 0);
    this->m_stats->iterations++;
    return this->working_set().size();
  }

  // FIXME move this to the template
  size_t end_iteration(std::vector<R>& parameters, std::vector<R>& action) override {
    auto result = end_iteration(wrap(parameters), wrap(action));
    return result;
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
  void set_max_size_qspace(int n) { m_max_size_qspace = n; }
  int get_max_size_qspace() const { return m_max_size_qspace; }
  void set_method(const std::string& method) { m_method = method; }
  std::string get_method() const { return m_method; }

  void set_options(const Options& options) override {
    SolverTemplate::set_options(options);
    auto opt = CastOptions::Optimize(options);
    if (opt.max_size_qspace)
      set_max_size_qspace(opt.max_size_qspace.value());
    if (opt.norm_thresh)
      set_norm_thresh(opt.norm_thresh.value());
    if (opt.svd_thresh)
      set_svd_thresh(opt.svd_thresh.value());
    if (opt.method)
      set_method(opt.method.value());
  }

  std::shared_ptr<Options> get_options() const override {
    auto opt = std::make_shared<OptimizeOptions>();
    opt->copy(*SolverTemplate::get_options());
    opt->max_size_qspace = get_max_size_qspace();
    opt->norm_thresh = get_norm_thresh();
    opt->svd_thresh = get_svd_thresh();
    opt->method = get_method();
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

  bool add_value(R& parameters, value_type value, R& residual) override {
    this->m_values.push(value);
    this->m_xspace->data[subspace::EqnData::value].resize({this->m_values.size(), 1});
    this->m_xspace->data[subspace::EqnData::value](this->m_values.size() - 1, 0) =
        value; // TODO find a less hacky way to inject value
    auto nwork = this->add_vector(parameters, residual);
    return nwork > 0; // TODO check this does the right thing
  }

  scalar_type value() const override { return this->m_values.top(); }

protected:
  // for non-linear problems, actions already contains the residual
  void construct_residual(const std::vector<int>& roots, const CVecRef<R>& params, const VecRef<R>& actions) override {}

  double m_norm_thresh = 1e-10; //!< vectors with norm less than threshold can be considered null.
  double m_svd_thresh = 1e-12;  //!< svd values smaller than this mark the null space
  int m_max_size_qspace = std::numeric_limits<int>::max(); //!< maximum size of Q space
  detail::DSpaceResetter<Q> m_dspace_resetter;             //!< resets D space
  bool m_hermiticity = true;                               //!< whether the problem is hermitian or not
  std::string m_method;                                    ///!< algorithm choice for solver
};

} // namespace molpro::linalg::itsolv
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_OPTIMIZE_H
