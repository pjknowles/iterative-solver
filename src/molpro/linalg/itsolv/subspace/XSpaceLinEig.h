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

template <class R, class Q, class P, typename ST>
class XSpaceLinEig : public XSpace<RSpace<R, Q, P>, QSpace<R, Q, P>, PSpace<R, P>, ST> {
  using XS = XSpace<RSpace<R, Q, P>, QSpace<R, Q, P>, PSpace<R, P>, ST>;

public:
  using typename XS::PS;
  using typename XS::QS;
  using typename XS::RS;
  using typename XS::scalar_type;
  using XS::data;
  XSpaceLinEig() = default;

  void check_conditioning(RS& rs, QS& qs, PS& ps) override {
    xspace::check_conditioning(*this, rs, qs, ps, m_svd_stability_threshold);
  }

  void solve(const IterativeSolver<R, Q, P>& solver) override {
    assert("XSpaceLinEig can only be used with LinearEigensystem solver");
  };

  void solve(const LinearEigensystem<R, Q, P>& solver) {
    auto& h = data[EqnData::H];
    auto& s = data[EqnData::S];
    if (m_hermitian)
      util::matrix_symmetrize(h);
    auto dim = h.rows();
    auto evec = std::vector<scalar_type>{};
    itsolv::eigenproblem(evec, m_eval, h.data(), s.data(), dim, m_hermitian, m_svd_solver_threshold, 0);
    auto n_solutions = evec.size() / dim;
    auto full_matrix = Matrix<scalar_type>{std::move(evec), {dim, n_solutions}};
    auto nroots = solver.n_roots();
    assert(solver.n_roots() == m_roots_in_subspace.size());
    assert(n_solutions >= solver.n_roots());
    m_eval.resize(nroots);
    m_evec.resize({dim, nroots});
    m_evec.slice() = full_matrix.slice({0, 0}, {dim, nroots});
    auto root_subspace = Matrix<double>({nroots, nroots});
    for (size_t i = 0; i < m_roots_in_subspace.size(); ++i)
      root_subspace.row(i) = m_evec.row(m_roots_in_subspace[i]);
    m_roots = util::eye_order(root_subspace);
  }

  const std::vector<scalar_type>& eigenvalues() const override { return m_eval; };

  const Matrix<scalar_type>& solutions() const override { return m_evec; };

  const std::vector<size_t>& roots() const override { return m_roots; };

  void build_subspace(RS& rs, QS& qs, PS& ps) override {
    m_dim = xspace::Dimensions(ps.size(), qs.size(), rs.size());
    xspace::build_subspace_H_S(data, rs.data, qs.data, qs.qr, qs.rq, ps.data, m_dim);
    // TODO make sure that there are checks to ensure converged and working set never overlap
    m_roots_in_subspace = xspace::roots_in_subspace(qs.converged_solutions(), rs.working_set(), m_dim.oQ, m_dim.oR);
  }

  const xspace::Dimensions& dimensions() const override { return m_dim; }

protected:
  xspace::Dimensions m_dim;
  bool m_hermitian = false; //!< whether the matrix is Hermitian
  double m_svd_stability_threshold =
      1.0e-4; //!< singular values of overlap matrix larger than this constitute a stable subspace
  std::map<size_t, std::vector<double>> m_solutions; //!< solutions mapped to root index
  double m_svd_solver_threshold = 1.0e-14;           //!< threshold to remove the null space during solution
  Matrix<scalar_type> m_evec;                        //!< eigenvectors stored as columns with ascending eigenvalue
  std::vector<scalar_type> m_eval;                   //!< eigenvalues in ascending order
  std::vector<size_t> m_roots;                       //!< for each eigenvector stores corresponding root index
  std::vector<size_t> m_roots_in_subspace; //!< indices of roots in the full subspace. Includes converged roots from
                                           //!< QSpace and working set from RSpace
};

} // namespace subspace
} // namespace itsolv
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_XSPACELINEIG_H
