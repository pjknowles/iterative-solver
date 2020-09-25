#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_LINEAREIGENSYSTEMA_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_LINEAREIGENSYSTEMA_H
#include <iterator>
#include <molpro/linalg/itsolv/IterativeSolverTemplate.h>
#include <molpro/linalg/itsolv/Logger.h>
#include <molpro/linalg/itsolv/subspace/XSpaceLinEig.h>

namespace molpro {
namespace linalg {
namespace itsolv {

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
class LinearEigensystemA : public IterativeSolverTemplate<
                               LinearEigensystem<R, Q, P>,
                               subspace::XSpaceLinEig<R, Q, P, typename LinearEigensystem<R, Q, P>::scalar_type>> {
public:
  using SolverTemplate =
      IterativeSolverTemplate<LinearEigensystem<R, Q, P>,
                              subspace::XSpaceLinEig<R, Q, P, typename LinearEigensystem<R, Q, P>::scalar_type>>;
  using typename SolverTemplate::scalar_type;

  explicit LinearEigensystemA(const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers,
                              std::shared_ptr<Logger> logger_ = std::make_shared<Logger>())
      : SolverTemplate(subspace::RSpace<R, Q, P>(handlers, logger_), subspace::QSpace<R, Q, P>(handlers, logger_),
                       subspace::PSpace<R, P>(), subspace::XSpaceLinEig<R, Q, P, scalar_type>(logger_), handlers,
                       std::make_shared<Statistics>(), logger_),
        logger(std::move(logger_)) {}

  void set_convergence_threshold(double threshold) { this->m_convergence_threshold = threshold; }

  std::vector<scalar_type> eigenvalues() const override { return this->m_xspace.eigenvalues(); };

  std::vector<scalar_type> working_set_eigenvalues() const {
    auto eval = std::vector<scalar_type>{};
    for (auto i : this->working_set()) {
      eval.emplace_back(this->m_xspace.eigenvalues().at(i));
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