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

namespace xspace {
struct Dimensions {
  Dimensions() : Dimensions(0, 0, 0) {}
  Dimensions(size_t np, size_t nq, size_t nr) : nP(np), nQ(nq), nR(nr), nX(nP + nQ + nR), oP(0), oQ(nP), oR(oQ + nQ) {}
  size_t nP, nQ, nR, nX, oP, oQ, oR;
};

//! Combines data from P, Q, and R subspaces to form the X subspace
void build_subspace_HS(SubspaceData& xx, const SubspaceData& rr, const SubspaceData& qq, const SubspaceData& qr,
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

void build_subspace_RHS(SubspaceData& xx, const SubspaceData& rr, const SubspaceData& qq, const SubspaceData& qr,
                        const SubspaceData& rq, const SubspaceData& pp, const Dimensions& d) {
  auto e = EqnData::rhs;
  // TODO copy data from subspaces
}

} // namespace xspace

//! Base class for XSpace solvers
template <class R, class Q, class P>
class XSpace {
public:
  using RS = RSpace<R, Q>;
  using QS = QSpace<R, Q, P>;
  using PS = PSpace<R, P>;

  SubspaceData data = null_data<EqnData::H, EqnData::S>();

  //! Make sure that the subspace is well conditioned, by merging or deleting elements of Q subspace
  void check_conditioning(RSpace<R, Q>& rs, QSpace<R, Q, P>& qs, PSpace<R, P>& ps) {
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
      build_subspace(rs, qs, ps);
    }
  };

  // Derived classes should implement a strategy for a particular solver
  // void solve(const IterativeSolver<R, Q, P>& solver);

  //! Access solution corresponding to root i
  virtual const std::vector<double>& solution(size_t i) const = 0;

  //! Number of vectors forming the subspace
  size_t size() { return data.at(EqnData::H).rows(); }

  void build_subspace(RSpace<R, Q>& rs, QSpace<R, Q, P>& qs, PSpace<R, P>& ps) {
    auto nP = ps.data.at(EqnData::H).rows();
    auto nQ = qs.data.at(EqnData::H).rows();
    auto nR = rs.data.at(EqnData::H).rows();
    m_dim = Dimensions(nP, nQ, nR);
    if (data.find(EqnData::rhs) == data.end())
      xspace::build_subspace_HS(data, rs.data, qs.data, qs.qr(), rs.rq(), ps.data, m_dim);
    else
      xspace::build_subspace_RHS(data, rs.data, qs.data, qs.qr(), rs.rq(), ps.data, m_dim);
  }

protected:
  xspace::Dimensions m_dim;
  double m_svd_stability_threshold =
      1.0e-4; //!< singular values of overlap matrix larger than this constitute a stable subspace
};

} // namespace subspace
} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_XSPACE_H