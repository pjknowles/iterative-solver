#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_PSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_PSPACE_H
#include <molpro/linalg/itsolv/subspace/SubspaceData.h>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace subspace {

template <class Rt, class Pt>
struct PSpace {
  using R = Rt;
  using P = Pt;
  using VecRefP = std::vector<std::reference_wrapper<P>>;

  SubspaceData data = null_data<EqnData::H, EqnData::S>(); //!< PxP block of subspace data
  VecRefP params() const { return {}; }
  VecRefP actions() const { return {}; }

  size_t size() const { return 0; }
};

} // namespace subspace
} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_PSPACE_H
