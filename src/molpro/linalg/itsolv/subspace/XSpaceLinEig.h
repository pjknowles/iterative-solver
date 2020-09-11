#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_XSPACELINEIG_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_XSPACELINEIG_H
#include <molpro/linalg/itsolv/helper.h>
#include <molpro/linalg/itsolv/subspace/XSpace.h>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace subspace {

template <class R, class Q, class P>
class XSpaceLinEig : public XSpace<R, Q, P> {
public:
  void solve(const LinearEigensystem<R, Q, P>& solver) {
    // itsolv::eigenproblem(m_evec_xx, m_eval_xx, m_h_xx, m_s_xx, m_n_x, m_hermitian, m_svdThreshold, m_verbosity);
  }

  //! Return solution vector for root i
  const std::vector<double>& solution(size_t i) const override { return m_solutions.at(i); };

protected:
  std::map<size_t, std::vector<double>> m_solutions; //!< solutions mapped to root index
  double m_svd_solver_threshold = 1.0e-14;           //!< threshold to remove the null space during solution
};

} // namespace subspace
} // namespace itsolv
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_XSPACELINEIG_H
