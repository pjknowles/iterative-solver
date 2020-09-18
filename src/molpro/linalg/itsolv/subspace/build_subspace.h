#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_BUILD_SUBSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_BUILD_SUBSPACE_H
#include <molpro/linalg/itsolv/subspace/XSpace.h>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace subspace {
namespace xspace {
//! Combines data from P, Q, and R subspaces to form the X subspace
void build_subspace_H_S(SubspaceData& xx, const SubspaceData& rr, const SubspaceData& qq, const SubspaceData& qr,
                        const SubspaceData& rq, const SubspaceData& pp, const Dimensions& d) {
  for (auto e : {EqnData::H, EqnData::S}) {
    xx.at(e).resize({d.nX, d.nX});
    xx.at(e).slice({d.oP, d.oP}, {d.nP, d.nP}) = pp.at(e);
    xx.at(e).slice({d.oQ, d.oQ}, {d.oQ + d.nQ, d.oQ + d.nQ}) = qq.at(e);
    xx.at(e).slice({d.oQ, d.oR}, {d.oQ + d.nQ, d.oR + d.nR}) = qr.at(e);
    xx.at(e).slice({d.oR, d.oQ}, {d.oR + d.nR, d.oQ + d.nQ}) = rq.at(e);
    xx.at(e).slice({d.oR, d.oR}, {d.oR + d.nR, d.oR + d.nR}) = rr.at(e);
  }
}

//! Combines data for the RHS equations in LinearEquations into the X subspace.
void build_subspace_RHS(SubspaceData& xx, const SubspaceData& rr, const SubspaceData& qq, const SubspaceData& qr,
                        const SubspaceData& rq, const SubspaceData& pp, const Dimensions& d) {
  auto e = EqnData::rhs;
  // TODO copy data from subspaces
}

} // namespace xspace
} // namespace subspace
} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_BUILD_SUBSPACE_H
