#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_QSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_QSPACE_H
#include <molpro/linalg/itsolv/detail/PSpace.h>
#include <molpro/linalg/itsolv/detail/RSpace.h>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace detail {

void update_qspace(SubspaceData<EqnData::H, EqnData::S>& qs, const SubspaceData<EqnData::H, EqnData::S>& rs) {}

template <class R, class Q, class P>
void update_qspace(SubspaceData<EqnData::H, EqnData::S, EqnData::rhs>& qs,
                   const SubspaceData<EqnData::H, EqnData::S, EqnData::rhs>& rs, LinearEigensystem<R, Q, P>& solver) {}

template <class R, class Q>
struct QSpace {
  SubspaceData<EqnData::H, EqnData::S> subspace;

  template <class P>
  void update(const RSpace<R>& rs, IterativeSolver<R, Q, P>& solver) {
    detail::update_qspace(subspace, rs.subspace);
  }
};

template <class R, class Q>
struct QSpaceLE {
  SubspaceData<EqnData::H, EqnData::S, EqnData::rhs> subspace;

  template <class P>
  void update(const RSpace<R>& rs, LinearEigensystem<R, Q, P>& solver) {
    detail::update_qspace(rs.subspace, solver);
  }
};

} // namespace detail
} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_QSPACE_H
