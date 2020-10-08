#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SUBSPACESOLVERLINEIG_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SUBSPACESOLVERLINEIG_H
#include <molpro/linalg/itsolv/subspace/Matrix.h>
#include <molpro/linalg/itsolv/subspace/XSpaceI.h>
#include <molpro/linalg/itsolv/subspace/check_conditioning.h>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace subspace {
template <typename VT, typename VTabs>
class SubspaceSolverLinEig {
public:
  using value_type = VT;
  using value_type_abs = VTabs;
  explicit SubspaceSolverLinEig(std::shared_ptr<Logger> logger) : m_logger(std::move(logger)) {}

  template <class R, class Q, class P>
  void check_conditioning(XSpaceI<R, Q, P>& xspace) {
    auto nX_on_entry = xspace.dimensions().nX;
    auto nX = xspace.dimensions().nX;
    m_logger->msg("SubspaceSolverLinEig::check_conditioning size of x space = " + std::to_string(nX), Logger::Trace);
    if (m_logger->data_dump) {
      m_logger->msg("on entry", Logger::Info);
      m_logger->msg("Sxx = " + as_string(xspace.data[EqnData::S]), Logger::Info);
      m_logger->msg("Hxx = " + as_string(xspace.data[EqnData::H]), Logger::Info);
    }
    xspace::check_conditioning_gram_schmidt(xspace, m_lin_trans, m_norm_stability_threshold, *m_logger);
    nX = xspace.dimensions().nX;
    m_logger->msg("size of x space after conditioning = " + std::to_string(nX), Logger::Debug);
    if (m_logger->data_dump && nX != nX_on_entry) {
      m_logger->msg("Sxx = " + as_string(xspace.data[EqnData::S]), Logger::Info);
      m_logger->msg("Hxx = " + as_string(xspace.data[EqnData::H]), Logger::Info);
    }
  }

  template <class R, class Q, class P>
  void solve(XSpaceI<R, Q, P>& xspace, const size_t nroots_max) {
    m_logger->msg("SubspaceSolverLinEig::solve", Logger::Trace);
    check_conditioning(xspace);
    // apply linear transformation
    auto& h = xspace.data[EqnData::H];
    auto& s = xspace.data[EqnData::S];
    auto dim = h.rows();
    auto evec = std::vector<value_type>{};
    itsolv::eigenproblem(evec, m_eigenvalues, h.data(), s.data(), dim, m_hermitian, m_svd_solver_threshold, 0);
    auto n_solutions = evec.size() / dim;
    auto full_matrix = Matrix<value_type>{std::move(evec), {n_solutions, dim}};
    auto nroots = std::min(nroots_max, n_solutions);
    m_eigenvalues.resize(nroots);
    m_solutions.resize({nroots, dim});
    m_solutions.slice() = full_matrix.slice({0, 0}, {nroots, dim});
    if (m_logger->data_dump) {
      m_logger->msg("eigenvalues = ", begin(m_eigenvalues), end(m_eigenvalues), Logger::Debug);
      m_logger->msg("eigenvectors = " + as_string(m_solutions), Logger::Info);
    }
  }

  const Matrix<value_type>& solutions() const { return m_solutions; }
  const std::vector<value_type>& eigenvalues() const { return m_eigenvalues; }

  size_t size() const { return m_solutions.rows(); }

protected:
  Matrix<value_type> m_solutions;        //!< solution matrix with row vectors
  std::vector<value_type> m_eigenvalues; //!< eigenvalues
  Matrix<value_type> m_lin_trans;        //!< linear transformation to a well conditioned subspace
  std::shared_ptr<Logger> m_logger;
  value_type_abs m_norm_stability_threshold; //!< norm threshold for Gram Schmidt orthogonalisation
  value_type_abs m_svd_solver_threshold;     //!< threshold to select null space during SVD in eigenproblem
  bool m_hermitian = true;                   //!< flags the matrix as Hermitian
};

} // namespace subspace
} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SUBSPACESOLVERLINEIG_H
