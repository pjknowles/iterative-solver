#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_XSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_XSPACE_H
#include <molpro/linalg/itsolv/detail/PSpace.h>
#include <molpro/linalg/itsolv/detail/QSpace.h>
#include <molpro/linalg/itsolv/detail/RSpace.h>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace detail {

//! Base class for XSpace solvers
template <class RSpace, class QSpace, class PSpace>
struct XSpace {
  container::Matrix subspace;
  //! Solver the subspace problem
  void solve() = 0;
  void check_conditioning(RSpace& rs, QSpace& qs, PSpace& ps) = 0;

  const container::Matrix& solutions() {}

  const std::vector<double>& solution(size_t i) {}

  size_t size() {}

  void build_subspace(RSpace& rs, QSpace& qs, PSpace& ps) {
    auto nP = ps.size();
    auto nQ = qs.size();
    auto nR = rs.size();
    size_t oP = 0;
    auto oQ = oP;
    auto oR = oQ + nQ;
    auto nX = nP + nQ + nR;
    subspace.resize(nX, nX);
    subspace.slice({oP, oP}, {nP, nP}) = ps.subspace;
    subspace.slice({oQ, oQ}, {oQ + nQ, oQ + nQ}) = qs.subspace();
    subspace.slice({oQ, oR}, {oQ + nQ, oR + nR}) = qs.qr();
    subspace.slice({oR, oQ}, {oR + nR, oQ + nQ}) = qs.rq();
    subspace.slice({oR, oR}, {oR + nR, oR + nR}) = rs.subspace();
  }
};

} // namespace detail
} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_XSPACE_H
