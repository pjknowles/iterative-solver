#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_PSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_PSPACE_H
#include <molpro/linalg/itsolv/subspace/SubspaceData.h>
#include <molpro/linalg/itsolv/wrap.h>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace subspace {

template <class Rt, class Pt>
struct PSpace {
  using R = Rt;
  using P = Pt;

  SubspaceData data = null_data<EqnData::H, EqnData::S>(); //!< PxP block of subspace data

  CVecRef<P> params() const { return {}; }
  CVecRef<P> cparams() const { return params(); }
  VecRef<P> params() { return {}; }
  CVecRef<P> actions() const { return {}; }
  CVecRef<P> cactions() const { return actions(); }
  VecRef<P> actions() { return {}; }

  size_t size() const { return 0; }
};

} // namespace subspace
} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_PSPACE_H
