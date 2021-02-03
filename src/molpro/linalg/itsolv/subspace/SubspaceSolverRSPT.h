#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SUBSPACESOLVERRSPT_H_
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SUBSPACESOLVERRSPT_H_
#include "SubspaceSolverLinEig.h"
namespace molpro::linalg::itsolv::subspace {
template <class RT, class QT, class PT>
class SubspaceSolverRSPT : public SubspaceSolverLinEig<RT, QT, PT> {
public:
  using value_type = typename ISubspaceSolver<RT, QT, PT>::value_type;
  using value_type_abs = typename ISubspaceSolver<RT, QT, PT>::value_type_abs;
  using R = typename ISubspaceSolver<RT, QT, PT>::R;
  using Q = typename ISubspaceSolver<RT, QT, PT>::Q;
  using P = typename ISubspaceSolver<RT, QT, PT>::P;

  explicit SubspaceSolverRSPT(std::shared_ptr<Logger> logger) : SubspaceSolverLinEig<RT, QT, PT>(std::move(logger)) {}

  void solve(IXSpace<R, Q, P>& xspace, const size_t nroots_max) override {
    SubspaceSolverLinEig<RT, QT, PT>::solve_eigenvalue(xspace, nroots_max);
    std::cout << "variational subspace solution " << as_string(this->m_solutions) << std::endl;
    std::cout << "rows " << this->m_solutions.rows() << std::endl;
    std::cout << "cols " << this->m_solutions.cols() << std::endl;
    auto n = this->m_solutions.cols();
    this->m_solutions.slice().fill(0);
    if (n > 1)
      this->m_solutions(0, n - 1) = 1;
    std::cout << "perturbational subspace solution " << as_string(this->m_solutions) << std::endl;
  }
};
} // namespace molpro::linalg::itsolv::subspace
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SUBSPACESOLVERRSPT_H_
