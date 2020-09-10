#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_XSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_XSPACE_H
#include <molpro/linalg/itsolv/subspace/PSpace.h>
#include <molpro/linalg/itsolv/subspace/QSpace.h>
#include <molpro/linalg/itsolv/subspace/RSpace.h>

#include <utility>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace subspace {

//! Combines data from P, Q, and R subspaces to form the X subspace
void build_subspace(SubspaceData& xs, const SubspaceData& rs, const SubspaceData& qq, const SubspaceData& qr,
                    const SubspaceData& rq, const SubspaceData& ps) {
  auto nP = ps.at(EqnData::H).rows();
  auto nQ = qq.at(EqnData::H).rows();
  auto nR = rs.at(EqnData::H).rows();
  size_t oP = 0;
  auto oQ = nP;
  auto oR = oQ + nQ;
  auto nX = nP + nQ + nR;
  for (auto d : {EqnData::H, EqnData::S}) {
    xs.at(d).resize({nX, nX});
    xs.at(d).slice({oP, oP}, {nP, nP}) = ps.at(d);
    xs.at(d).slice({oQ, oQ}, {oQ + nQ, oQ + nQ}) = qq.at(d);
    xs.at(d).slice({oQ, oR}, {oQ + nQ, oR + nR}) = qr.at(d);
    xs.at(d).slice({oR, oQ}, {oR + nR, oQ + nQ}) = rq.at(d);
    xs.at(d).slice({oR, oR}, {oR + nR, oR + nR}) = rs.at(d);
  }
}

//! Base class for XSpace solvers
template <class Rs, class Qs, class Ps>
struct XSpace {
  using RS = Rs;
  using QS = Qs;
  using PS = Ps;
  SubspaceData data = null_data<EqnData::H, EqnData::S>();

  //! Make sure that the subspace is well conditioned. Derived classes should implement a strategy
  virtual void check_conditioning(RS& rs, QS& qs, PS& ps) = 0;
  //! Solve the subspace problem
  virtual void solve() = 0;

  const Matrix<double>& solutions() {}

  std::vector<double> solution(size_t i) {}

  //! Number of vectors forming the subspace
  size_t size() { return data.at(EqnData::H).rows(); }

  void build_subspace(const RS& rs, const QS& qs, const PS& ps) {
    subspace::build_subspace(data, rs.data, qs.data, qs.qr(), qs.rq(), ps.data);
  }
};

} // namespace detail
} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_XSPACE_H