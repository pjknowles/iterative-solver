#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_OPTIMIZEBFGS_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_OPTIMIZEBFGS_H
#include <molpro/linalg/itsolv/CastOptions.h>
#include <molpro/linalg/itsolv/DSpaceResetter.h>
#include <molpro/linalg/itsolv/IterativeSolverTemplate.h>
#include <molpro/linalg/itsolv/propose_rspace.h>
#include <molpro/linalg/itsolv/subspace/SubspaceSolverOptBFGS.h>
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
template <class R, class Q, class P = std::map<size_t, typename R::value_type>>
class OptimizeBFGS : public IterativeSolverTemplate<Optimize, R, Q, P> {
public:
  using SolverTemplate = IterativeSolverTemplate<Optimize, R, Q, P>;
  using SolverTemplate ::report;
  using typename SolverTemplate::scalar_type;
  using typename SolverTemplate::value_type;
  using typename SolverTemplate::value_type_abs;
  using SubspaceSolver = subspace::SubspaceSolverOptBFGS<R, Q, P>;

  explicit OptimizeBFGS(const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers,
                        const std::shared_ptr<Logger>& logger_ = std::make_shared<Logger>())
      : SolverTemplate(std::make_shared<subspace::XSpace<R, Q, P>>(handlers, logger_),
                       std::static_pointer_cast<subspace::ISubspaceSolver<R, Q, P>>(
                           std::make_shared<subspace::SubspaceSolverOptBFGS<R, Q, P>>(logger_)),
                       handlers, std::make_shared<Statistics>(), logger_),
        logger(logger_) {}

  size_t end_iteration(const VecRef<R>& parameters, const VecRef<R>& action) override {
    this->solution_params(this->m_working_set, parameters);
    if (this->m_errors.front() < this->m_convergence_threshold) {
      this->m_working_set.clear();
      return 0;
    }
    this->m_working_set.assign(1, 0);
    auto& z = action.front();
    const auto& xspace = this->m_xspace;
    auto& xdata = xspace->data;
    const auto& H = xdata[subspace::EqnData::H];
    const auto& S = xdata[subspace::EqnData::S];
    const auto& q = xspace->paramsq();
    const auto& u = xspace->actionsq();

    auto signal = solver_signal();
    //    std::cout << "end_iteration signal=" << signal << std::endl;
    if (signal == 0) { // action is expected to hold the preconditioned residual
      for (int a = m_alpha.size() - 1; a >= 0; a--) {
        auto beta = -(this->m_handlers->qr().dot(z, u[a]) - this->m_handlers->qr().dot(z, u[a + 1])) /
                    (H(a, a) - H(a, a + 1) - H(a + 1, a) + H(a + 1, a + 1));
        //        std::cout << "alpha[" << a << "] = " << m_alpha[a] << std::endl;
        //        std::cout << "beta[" << a << "] = " << beta << std::endl;
        this->m_handlers->qr().axpy(-m_alpha[a] + beta, q[a], z);
        this->m_handlers->qr().axpy(+m_alpha[a] - beta, q[a + 1], z);
      }
      this->m_handlers->rr().axpy(1, z, parameters.front());
    } else if (signal == 1) { // residual not used, simply leave parameters alone
    }
    this->m_stats->iterations++;
    return this->errors().front() < this->m_convergence_threshold ? 0 : 1;
  }

  size_t end_iteration(std::vector<R>& parameters, std::vector<R>& action) override {
    auto result = end_iteration(wrap(parameters), wrap(action));
    return result;
  }

  //! Set a limit on the maximum size of Q space. This does not include the size of the working space (R) and the D
  //! space
  void set_max_size_qspace(int n) { m_max_size_qspace = n; }
  int get_max_size_qspace() const { return m_max_size_qspace; }

  void set_options(const Options& options) override {
    SolverTemplate::set_options(options);
    if (auto opt2 = dynamic_cast<const molpro::linalg::itsolv::OptimizeBFGSOptions*>(&options)) {
      auto opt = CastOptions::OptimizeBFGS(options);
      if (opt.max_size_qspace)
        set_max_size_qspace(opt.max_size_qspace.value());
      if (opt.strong_Wolfe)
        m_strong_Wolfe = opt.strong_Wolfe.value();
      if (opt.Wolfe_1)
        m_Wolfe_1 = opt.Wolfe_1.value();
      if (opt.Wolfe_2)
        m_Wolfe_2 = opt.Wolfe_2.value();
      if (opt.linesearch_tolerance)
        m_linesearch_tolerance = opt.linesearch_tolerance.value();
      if (opt.linesearch_grow_factor)
        m_linesearch_grow_factor = opt.linesearch_grow_factor.value();
    }
    if (auto opt2 = dynamic_cast<const molpro::linalg::itsolv::OptimizeSDOptions*>(&options)) {
      auto opt = CastOptions::OptimizeSD(options);
    }
  }

  std::shared_ptr<Options> get_options() const override {
    auto opt = std::make_shared<OptimizeBFGSOptions>();
    opt->copy(*SolverTemplate::get_options());
    opt->max_size_qspace = get_max_size_qspace();
    opt->strong_Wolfe = m_strong_Wolfe;
    opt->Wolfe_1 = m_strong_Wolfe;
    opt->Wolfe_2 = m_Wolfe_2;
    opt->linesearch_tolerance = m_linesearch_tolerance;
    opt->linesearch_grow_factor = m_linesearch_grow_factor;
    return opt;
  }

  void report(std::ostream& cout) const override {
    SolverTemplate::report(cout);
    cout << "value " << value() << ", errors " << std::scientific;
    auto& err = this->m_errors;
    std::copy(begin(err), end(err), std::ostream_iterator<value_type_abs>(molpro::cout, ", "));
    cout << std::defaultfloat << std::endl;
  }
  std::shared_ptr<Logger> logger;

protected:
  std::vector<double> m_alpha;

public:
  bool add_value(R& parameters, value_type value, R& residual) override {
    using namespace subspace;
    auto& xspace = this->m_xspace;
    auto& xdata = xspace->data;
    const auto& H = xdata[EqnData::H];
    const auto& S = xdata[EqnData::S];
    auto& Value = xdata[EqnData::value];
    const int n = this->m_xspace->dimensions().nX;
    auto oldValue = Value;
    Value.resize({n + 1, 1});
    if (n > 0)
      Value.slice({1, 0}, {n + 1, 1}) = oldValue.slice();
    Value(0, 0) = value;
    //    std::cout << "updated value to" << as_string(Value) << std::endl;
    auto nwork = this->add_vector(parameters, residual);
    if (solver_signal() == 1) // line search required
      return false;
    m_alpha.resize(n);
    const auto& q = xspace->paramsq();
    const auto& u = xspace->actionsq();
    //    this->m_errors.front() = std::sqrt(this->m_handlers->rr().dot(residual,residual));
    for (int a = 0; a < n; a++) {
      m_alpha[a] = (this->m_handlers->qr().dot(residual, q[a]) - this->m_handlers->qr().dot(residual, q[a + 1])) /
                   (H(a, a) - H(a, a + 1) - H(a + 1, a) + H(a + 1, a + 1));
      //      std::cout << "alpha[" << a << "] = " << m_alpha[a] << std::endl;
      this->m_handlers->qr().axpy(-m_alpha[a], u[a], residual);
      this->m_handlers->qr().axpy(m_alpha[a], u[a + 1], residual);
    }
    return nwork > 0;
  }

  size_t end_iteration(R& parameters, R& actions) override {
    auto wparams = std::vector<std::reference_wrapper<R>>{std::ref(parameters)};
    auto wactions = std::vector<std::reference_wrapper<R>>{std::ref(actions)};
    return end_iteration(wparams, wactions);
  }

  scalar_type value() const override { return this->m_xspace->data[subspace::EqnData::value](0, 0); }

protected:
  int solver_signal() const {
    int signal = 0;
    using namespace subspace;
    auto& xspace = this->m_xspace;
    auto& xdata = xspace->data;
    if (xdata.find(EqnData::signals) != xdata.end() && xdata[EqnData::signals].empty())
      signal = xdata[EqnData::signals](0, 0);
    //    molpro::cout << "signal " << signal << std::endl;
    return signal;
  }
  // for non-linear problems, actions already contains the residual
  void construct_residual(const std::vector<int>& roots, const CVecRef<R>& params, const VecRef<R>& actions) override {}

  int m_max_size_qspace = std::numeric_limits<int>::max(); //!< maximum size of Q space
  bool m_strong_Wolfe = true;                              //!< Whether to use strong or weak Wolfe conditions
  double m_Wolfe_1 = 1e-4; //!< Acceptance parameter for function value; recommended value Nocedal and Wright p142
  double m_Wolfe_2 = 0.9;  //!< Acceptance parameter for function gradient; recommended value Nocedal and Wright p142
  double m_linesearch_tolerance = .2; //!< If the predicted line search is within tolerance of 1, don't bother taking it
  double m_linesearch_grow_factor =
      3; //!< If the predicted line search step is extrapolation, limit the step to this factor times the current step
};

} // namespace molpro::linalg::itsolv
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_OPTIMIZEBFGS_H
