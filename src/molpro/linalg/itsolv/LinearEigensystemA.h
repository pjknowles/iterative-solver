#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_LINEAREIGENSYSTEMA_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_LINEAREIGENSYSTEMA_H
#include <molpro/linalg/itsolv/IterativeSolverTemplate.h>
#include <molpro/linalg/itsolv/detail/XSpace.h>

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
class LinearEigensystemA
    : public IterativeSolverTemplate<LinearEigensystem<R, Q, P>,
                                     detail::XSpace<detail::RSpace<R>, detail::QSpace<R, Q>, detail::PSpace<R, P>>> {
public:
  using typename LinearEigensystem<R, Q, P>::scalar_type;

  std::vector<scalar_type> eigenvalues() const override { return this->m_xspace.eigenvalues(); };
};

} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_LINEAREIGENSYSTEMA_H