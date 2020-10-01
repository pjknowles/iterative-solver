#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_XSPACELINEIG_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_XSPACELINEIG_H
#include <cassert>
#include <molpro/linalg/itsolv/helper.h>
#include <molpro/linalg/itsolv/subspace/CSpace.h>
#include <molpro/linalg/itsolv/subspace/XSpace.h>
#include <molpro/linalg/itsolv/subspace/build_subspace.h>
#include <molpro/linalg/itsolv/subspace/check_conditioning.h>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace subspace {

template <class R, class Q, class P, typename ST>
class XSpaceLinEig : public XSpace<RSpace<R, Q, P>, QSpace<R, Q, P>, PSpace<R, P>, CSpace<R, Q, P, ST>, ST> {
  using XS = XSpace<RSpace<R, Q, P>, QSpace<R, Q, P>, PSpace<R, P>, CSpace<R, Q, P, ST>, ST>;

public:
  using typename XS::CS;
  using typename XS::PS;
  using typename XS::QS;
  using typename XS::RS;
  using typename XS::scalar_type;
  using XS::data;
  explicit XSpaceLinEig(std::shared_ptr<Logger> logger) : m_logger(std::move(logger)){};

  void check_conditioning(RS& rs, QS& qs, PS& ps, CS& cs) override {
    auto nX_on_entry = m_dim.nX;
    m_logger->msg("XSpaceLinEig::check_conditioning size of x space = " + std::to_string(m_dim.nX), Logger::Trace);
    m_logger->msg("size of x space before conditioning = " + std::to_string(m_dim.nX), Logger::Debug);
    if (m_logger->data_dump) {
      m_logger->msg("on entry", Logger::Info);
      m_logger->msg("Sxx = " + as_string(data[EqnData::S]), Logger::Info);
      m_logger->msg("Hxx = " + as_string(data[EqnData::H]), Logger::Info);
    }
    xspace::check_conditioning_gram_schmidt(*this, rs, qs, ps, cs, m_subspace_transformation,
                                            m_norm_stability_threshold, *m_logger);
    m_logger->msg("size of x space after conditioning = " + std::to_string(m_dim.nX), Logger::Debug);
    if (m_logger->data_dump && m_dim.nX != nX_on_entry) {
      m_logger->msg("Sxx = " + as_string(data[EqnData::S]), Logger::Info);
      m_logger->msg("Hxx = " + as_string(data[EqnData::H]), Logger::Info);
    }
  }

  void solve(const IterativeSolver<R, Q, P>& solver) override {
    assert("XSpaceLinEig can only be used with LinearEigensystem solver");
  };

  void solve(const LinearEigensystem<R, Q, P>& solver) {
    m_logger->msg("XSpaceLinEig::solve", Logger::Trace);
    auto& h = data[EqnData::H];
    auto& s = data[EqnData::S];
    if (m_hermitian)
      util::matrix_symmetrize(h);
    auto dim = h.rows();
    auto evec = std::vector<scalar_type>{};
    itsolv::eigenproblem(evec, m_eval, h.data(), s.data(), dim, m_hermitian, m_svd_solver_threshold, 0);
    auto n_solutions = evec.size() / dim;
    auto full_matrix = Matrix<scalar_type>{std::move(evec), {n_solutions, dim}};
    auto nroots = std::min(solver.n_roots(), n_solutions);
    m_eval.resize(nroots);
    m_evec.resize({nroots, dim});
    m_evec.slice() = full_matrix.slice({0, 0}, {nroots, dim});
    if (m_logger->data_dump) {
      m_logger->msg("eigenvalues = ", begin(m_eval), end(m_eval), Logger::Debug);
      m_logger->msg("eigenvectors = " + as_string(m_evec), Logger::Info);
    }
  }

  const std::vector<scalar_type>& eigenvalues() const override { return m_eval; };

  const Matrix<scalar_type>& solutions() const override { return m_evec; };

  void build_subspace(RS& rs, QS& qs, PS& ps, CS& cs) override {
    m_dim = xspace::Dimensions(ps.size(), qs.size(), rs.size(), cs.size());
    xspace::build_subspace_H_S(data, ps.data, rs.data, qs.data, cs.data, qs.qr, qs.qc, qs.rq, qs.cq, m_dim);
  }

  const xspace::Dimensions& dimensions() const override { return m_dim; }

protected:
  std::shared_ptr<Logger> m_logger;
  xspace::Dimensions m_dim;
  bool m_hermitian = false; //!< whether the matrix is Hermitian
  double m_norm_stability_threshold =
      1.0e-5; //!< norm subspace vector after orhtogonalisation must be less than this to trigger removal
  double m_svd_solver_threshold = 1.0e-14; //!< threshold to remove the null space during solution
  Matrix<scalar_type>
      m_subspace_transformation;   //!< linear transformation of subspace vectors that leads to stable overlap
  Matrix<scalar_type> m_evec;      //!< eigenvectors stored as columns with ascending eigenvalue
  std::vector<scalar_type> m_eval; //!< eigenvalues in ascending order
};

} // namespace subspace
} // namespace itsolv
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_XSPACELINEIG_H
