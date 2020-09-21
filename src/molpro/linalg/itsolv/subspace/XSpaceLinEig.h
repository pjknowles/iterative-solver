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

  void check_conditioning(RS& rs, QS& qs, PS& ps) override {
    xspace::check_conditioning(*this, rs, qs, ps, m_svd_stability_threshold);
  }

  void solve(const IterativeSolver<R, Q, P>& solver) override {
    assert("XSpaceLinEig can only be used with LinearEigensystem solver");
  };

  void solve(const LinearEigensystem<R, Q, P>& solver) {
    auto& h = data[EqnData::H].data();
    auto& s = data[EqnData::H].data();
    if (m_hermitian)
      util::matrix_symmetrize(h);
    auto dim = h.rows();
    auto evec = std::vector<scalar_type>{};
    itsolv::eigenproblem(evec, m_eval, h, s, dim, m_hermitian, m_svd_solver_threshold, 0);
    auto n_solutions = evec.size() / dim;
    auto full_matrix = Matrix<scalar_type>{std::move(evec), {dim, n_solutions}};
    auto nR = solver.working_set().size();
    full_matrix.resize({dim, nR});
    m_evec = std::move(full_matrix);
    m_roots = util::eye_order(m_evec.slice({dim - nR, 0}, {dim, nR}));
  }

  const auto& eigenvalues() const { return m_eval; };

  const Matrix<scalar_type>& solution() const override { return m_evec; };

  const std::vector<size_t>& roots() const override { return m_roots; };

  void build_subspace(RS& rs, QS& qs, PS& ps) override {
    m_dim = Dimensions(ps.size(), qs.size(), rs.size());
    xspace::build_subspace_H_S(data, rs.data, qs.data, qs.qr(), rs.rq(), ps.data, m_dim);
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
  std::vector<size_t> m_roots;                       //!< for each solution stores corresponding root index
};

} // namespace subspace
} // namespace itsolv
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_XSPACELINEIG_H
