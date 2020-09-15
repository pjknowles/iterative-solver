#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_XSPACELINEIG_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_XSPACELINEIG_H
#include <cassert>
#include <molpro/linalg/itsolv/helper.h>
#include <molpro/linalg/itsolv/subspace/XSpace.h>
#include <molpro/linalg/itsolv/subspace/build_subspace.h>
#include <molpro/linalg/itsolv/subspace/check_conditioning.h>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace subspace {

template <class R, class Q, class P, typename scalar_type>
class XSpaceLinEig : public XSpace<RSpace<R, Q, P>, QSpace<R, Q, P>, PSpace<R, P>> {
  using XS = XSpace<RSpace<R, Q, P>, QSpace<R, Q, P>, PSpace<R, P>>;

public:
  using typename XS::PS;
  using typename XS::QS;
  using typename XS::RS;
  using XS::data;

  void check_conditioning(RS& rs, QS& qs, PS& ps) override { xspace::check_conditioning(*this, rs, qs, ps); }

  void solve(const IterativeSolver<R, Q, P>& solver) override {
    assert("XSpaceLinEig can only be used with LinearEigensystem solver");
  };

  void solve(const LinearEigensystem<R, Q, P>& solver) {
    // itsolv::eigenproblem(m_evec_xx, m_eval_xx, m_h_xx, m_s_xx, m_n_x, m_hermitian, m_svdThreshold, m_verbosity);
  }

  std::vector<scalar_type> eigenvalues() const { return {}; };

  //! Return solution vector for root i
  const std::vector<double>& solution(size_t i) const override { return m_solutions.at(i); };

  void build_subspace(RS& rs, QS& qs, PS& ps) override {
    auto nP = ps.data.at(EqnData::H).rows();
    auto nQ = qs.data.at(EqnData::H).rows();
    auto nR = rs.data.at(EqnData::H).rows();
    m_dim = Dimensions(nP, nQ, nR);
    xspace::build_subspace_HS(data, rs.data, qs.data, qs.qr(), rs.rq(), ps.data, m_dim);
  }

protected:
  xspace::Dimensions m_dim;
  double m_svd_stability_threshold =
      1.0e-4; //!< singular values of overlap matrix larger than this constitute a stable subspace
  std::map<size_t, std::vector<double>> m_solutions; //!< solutions mapped to root index
  double m_svd_solver_threshold = 1.0e-14;           //!< threshold to remove the null space during solution
};

} // namespace subspace
} // namespace itsolv
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_XSPACELINEIG_H
