#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_QSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_QSPACE_H
#include <molpro/linalg/itsolv/detail/PSpace.h>
#include <molpro/linalg/itsolv/detail/RSpace.h>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace detail {

template <class R, class Q, class P>
void update_qspace(SubspaceData& qs, const SubspaceData& rs, IterativeSolver<R, Q, P>& solver) {}

template <class R, class Q, class P>
void update_qspace(SubspaceData& qs, const SubspaceData& rs, LinearEigensystem<R, Q, P>& solver) {}

template <class R, class Q>
struct QSpace {
  SubspaceData subspace = null_data<EqnData::H, EqnData::S>();

  template <class P>
  void update(const RSpace<R>& rs, IterativeSolver<R, Q, P>& solver) {
    detail::update_qspace(subspace, rs.subspace);
  }
};

template <class R, class Q>
struct QSpaceLE : public QSpace<R, Q> {
  using QSpace<R, Q>::subspace;
  QSpaceLE() : QSpace<R, Q>() { subspace = null_data<EqnData::H, EqnData::S>; }

  template <class P>
  void update(const RSpace<R>& rs, LinearEigensystem<R, Q, P>& solver) {
    QSpaceLE<R, Q>::update(rs, solver);
  }
};

} // namespace detail
} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_QSPACE_H
