#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_LINEAREIGENSYSTEMA_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_LINEAREIGENSYSTEMA_H
#include <molpro/linalg/itsolv/IterativeSolverTemplate.h>
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

  explicit LinearEigensystemA(std::shared_ptr<ArrayHandlers<R, Q, P>> handlers)
      : SolverTemplate(subspace::RSpace<R, Q, P>(handlers), subspace::QSpace<R, Q, P>(handlers),
                       subspace::PSpace<R, P>(), subspace::XSpaceLinEig<R, Q, P, scalar_type>(), std::move(handlers),
                       std ::make_shared<Statistics>()) {}

  void set_convergence_threshold(double threshold) { this->m_convergence_threshold = threshold; }

  std::vector<scalar_type> eigenvalues() const override { return this->m_xspace.eigenvalues(); };

  std::vector<scalar_type> working_set_eigenvalues() const {
    auto eval = std::vector<scalar_type>{};
    for (auto i : this->working_set()) {
      eval.emplace_back(this->m_xspace.eigenvalues().at(i));
    }
    return eval;
  }
};

} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_LINEAREIGENSYSTEMA_H