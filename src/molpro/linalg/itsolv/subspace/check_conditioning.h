#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_CHECK_CONDITIONING_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_CHECK_CONDITIONING_H
#include <molpro/linalg/itsolv/subspace/PSpace.h>
#include <molpro/linalg/itsolv/subspace/QSpace.h>
#include <molpro/linalg/itsolv/subspace/RSpace.h>
#include <molpro/linalg/itsolv/subspace/XSpace.h>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace subspace {
namespace xspace {

template <class R, class P, class Q>
void check_conditioning(XSpace<RSpace<R, Q, P>, QSpace<R, Q, P>, PSpace<R, P>>& xs, RSpace<R, Q, P>& rs,
                        QSpace<R, Q, P>& qs, PSpace<R, P>& ps) {
  // TODO implement and move to a free function
  // get list of candidates for modification
  bool stable = false;
  std::vector<size_t> candidates;
  while (!stable && !candidates.empty()) {
    // Form a list of all adjacent pairs that can be merged
    // do SVD, return the lowest singular value
    // if singular value is greater than threshold, than it is stable, break
    // form a vector of norms for all pairs
    // pass the pair with largest norm for merger to q space
    xs.build_subspace(rs, qs, ps);
  }
}

} // namespace xspace
} // namespace subspace
} // namespace itsolv
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_CHECK_CONDITIONING_H
