#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_LINEAREIGENSYSTEMA_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_LINEAREIGENSYSTEMA_H
#include <iterator>
#include <molpro/linalg/itsolv/IterativeSolverTemplate.h>
#include <molpro/linalg/itsolv/Logger.h>
#include <molpro/linalg/itsolv/propose_rspace.h>
#include <molpro/linalg/itsolv/subspace/SubspaceSolverLinEig.h>
#include <molpro/linalg/itsolv/subspace/XSpace.h>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace detail {
template <class R>
using SubspaceSolver = subspace::SubspaceSolverLinEig<typename array::ArrayHandler<R, R>::value_type,
                                                      typename array::ArrayHandler<R, R>::value_type_abs>;
}

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
class LinearEigensystemA
    : public IterativeSolverTemplate<LinearEigensystem<R, Q, P>, subspace::XSpace<R, Q, P>, detail::SubspaceSolver<R>> {
public:
  using SolverTemplate =
      IterativeSolverTemplate<LinearEigensystem<R, Q, P>, subspace::XSpace<R, Q, P>, detail::SubspaceSolver<R>>;
  using typename SolverTemplate::scalar_type;

  explicit LinearEigensystemA(const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers,
                              const std::shared_ptr<Logger>& logger_ = std::make_shared<Logger>())
      : SolverTemplate(subspace::XSpace<R, Q, P>(handlers, logger_), detail::SubspaceSolver<R>{logger}, handlers,
                       std::make_shared<Statistics>(), logger_),
        logger(logger_) {}

  /*!
   * \brief Proposes new parameters for the subspace from the preconditioned residuals.
   *
   * After add_vector solves the subspace problem it returns the new solution and residual. The residual should be
   * preconditioned, i.e. using Davidson method, to accelerate convergence. The updated residual is used in this
   * function to propose new Q space parameters orthonormal to the old space. They are returned in parameters so that
   * corresponding actions can be calculated and used in add_vector in the next iteration.
   *
   * @param parameters output new parameters for the subspace.
   * @param residual preconditioned residuals.
   * @return number of significant parameters to calculate the action for
   */
  size_t end_iteration(std::vector<R>& parameters, std::vector<R>& action) override {
    auto r_norm_thresh = 1.0e-14;
    //    this->m_working_set = detail::propose_rspace(*static_cast<LinearEigensystem<R, Q, P>*>(this), parameters,
    //    action,
    //                                                 this->m_xspace.pspace, this->m_xspace.qspace,
    //                                                 this->m_xspace.cspace, this->m_xspace, *this->m_handlers,
    //                                                 *this->m_logger, r_norm_thresh);
    this->m_stats->iterations++;
    return this->working_set().size();
  }

  //! Applies the Davidson preconditioner
  void precondition(std::vector<R>& parameters, std::vector<R>& action) const {}

  void set_convergence_threshold(double threshold) { this->m_convergence_threshold = threshold; }

  std::vector<scalar_type> eigenvalues() const override { return this->m_subspace_solver.eigenvalues(); }

  std::vector<scalar_type> working_set_eigenvalues() const {
    auto eval = std::vector<scalar_type>{};
    for (auto i : this->working_set()) {
      eval.emplace_back(this->m_subspace_solver.eigenvalues().at(i));
    }
    return eval;
  }

  void report() const override {
    SolverTemplate::report();
    molpro::cout << "eigenvalues ";
    auto ev = eigenvalues();
    molpro::cout << std::fixed << std::setprecision(14);
    std::copy(begin(ev), end(ev), std::ostream_iterator<scalar_type>(molpro::cout, ", "));
    molpro::cout << std::defaultfloat << std::endl;
  }

  std::shared_ptr<Logger> logger;
};

} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_LINEAREIGENSYSTEMA_H