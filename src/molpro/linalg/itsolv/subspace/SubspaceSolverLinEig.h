#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SUBSPACESOLVERLINEIG_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SUBSPACESOLVERLINEIG_H
#include <molpro/linalg/itsolv/subspace/Matrix.h>
#include <molpro/linalg/itsolv/subspace/SubspaceSolverI.h>
#include <molpro/linalg/itsolv/subspace/XSpaceI.h>

namespace molpro::linalg::itsolv::subspace {

template <class RT, class QT, class PT>
class SubspaceSolverLinEig : public SubspaceSolverI<RT, QT, PT> {
public:
  using value_type = typename SubspaceSolverI<RT, QT, PT>::value_type;
  using value_type_abs = typename SubspaceSolverI<RT, QT, PT>::value_type_abs;
  using R = typename SubspaceSolverI<RT, QT, PT>::R;
  using Q = typename SubspaceSolverI<RT, QT, PT>::Q;
  using P = typename SubspaceSolverI<RT, QT, PT>::P;

  explicit SubspaceSolverLinEig(std::shared_ptr<Logger> logger) : m_logger(std::move(logger)) {}

  void solve(XSpaceI<R, Q, P>& xspace, const size_t nroots_max) override {
    m_logger->msg("SubspaceSolverLinEig::solve", Logger::Trace);
    auto h = xspace.data[EqnData::H];
    auto s = xspace.data[EqnData::S];
    if (m_logger->data_dump) {
      m_logger->msg("S = " + as_string(s), Logger::Info);
      m_logger->msg("H = " + as_string(h, 15), Logger::Info);
    }
//    std::cout << "h "<<as_string(h)<<std::endl;
//    std::cout << "s "<<as_string(s)<<std::endl;
    auto dim = h.rows();
    auto evec = std::vector<value_type>{};
    int verbosity = m_logger->max_trace_level == Logger::Info ? 3 : 0;
    itsolv::eigenproblem(evec, m_eigenvalues, h.data(), s.data(), dim, m_hermitian, m_svd_solver_threshold, verbosity);
    auto n_solutions = evec.size() / dim;
    auto full_matrix = Matrix<value_type>{std::move(evec), {n_solutions, dim}};
    auto nroots = std::min(nroots_max, n_solutions);
    m_eigenvalues.resize(nroots);
    m_solutions.resize({nroots, dim});
    m_solutions.slice() = full_matrix.slice({0, 0}, {nroots, dim});
    m_errors.assign(size(), std::numeric_limits<value_type_abs>::max());
    if (m_logger->data_dump) {
      m_logger->msg("eigenvalues = ", begin(m_eigenvalues), end(m_eigenvalues), Logger::Debug, 10);
      m_logger->msg("eigenvectors = " + as_string(m_solutions), Logger::Info);
    }
  }

  //! Set error value for solution *root*
  void set_error(int root, value_type_abs error) override { m_errors.at(root) = error; }
  void set_error(const std::vector<int>& roots, const std::vector<value_type_abs>& errors) override {
    for (size_t i = 0; i < roots.size(); ++i)
      set_error(roots[i], errors[i]);
  }

  const Matrix<value_type>& solutions() const override { return m_solutions; }
  const std::vector<value_type>& eigenvalues() const override { return m_eigenvalues; }
  const std::vector<value_type_abs>& errors() const override { return m_errors; }

  //! Number of solutions
  size_t size() const override { return m_solutions.rows(); }

protected:
  Matrix<value_type> m_solutions;        //!< solution matrix with row vectors
  std::vector<value_type> m_eigenvalues; //!< eigenvalues
  std::vector<value_type_abs> m_errors;  //!< errors in subspace solutions
  std::shared_ptr<Logger> m_logger;

public:
  value_type_abs m_svd_solver_threshold = 1.0e-14; //!< threshold to select null space during SVD in eigenproblem
  bool m_hermitian = true;                         //!< flags the matrix as Hermitian
};

} // namespace molpro::linalg::itsolv::subspace

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SUBSPACESOLVERLINEIG_H
