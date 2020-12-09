#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SUBSPACESOLVERLINEIG_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SUBSPACESOLVERLINEIG_H
#include <molpro/linalg/itsolv/subspace/ISubspaceSolver.h>
#include <molpro/linalg/itsolv/subspace/IXSpace.h>
#include <molpro/linalg/itsolv/subspace/Matrix.h>

namespace molpro::linalg::itsolv::subspace {

/*!
 * @brief Solves subspace problem for linear eigenvalues and system of linear equations
 */
template <class RT, class QT, class PT>
class SubspaceSolverLinEig : public ISubspaceSolver<RT, QT, PT> {
public:
  using value_type = typename ISubspaceSolver<RT, QT, PT>::value_type;
  using value_type_abs = typename ISubspaceSolver<RT, QT, PT>::value_type_abs;
  using R = typename ISubspaceSolver<RT, QT, PT>::R;
  using Q = typename ISubspaceSolver<RT, QT, PT>::Q;
  using P = typename ISubspaceSolver<RT, QT, PT>::P;

  explicit SubspaceSolverLinEig(std::shared_ptr<Logger> logger) : m_logger(std::move(logger)) {}

  void solve(IXSpace<R, Q, P>& xspace, const size_t nroots_max) override {
    m_logger->msg("SubspaceSolverLinEig::solve", Logger::Trace);
    if (xspace.data.find(EqnData::rhs) == xspace.data.end() || xspace.data[EqnData::rhs].empty()) {
      solve_eigenvalue(xspace, nroots_max);
    } else {
      solve_linear_equations(xspace);
    }
  }

protected:
  void solve_eigenvalue(IXSpace<R, Q, P>& xspace, const size_t nroots_max) {
    m_logger->msg("SubspaceSolverLinEig::solve_eigenvalue", Logger::Trace);
    auto h = xspace.data[EqnData::H];
    auto s = xspace.data[EqnData::S];
    if (m_logger->data_dump) {
      m_logger->msg("S = " + as_string(s), Logger::Info);
      m_logger->msg("H = " + as_string(h, 15), Logger::Info);
      m_logger->msg("rhs = " + as_string(h, 15), Logger::Info);
    }
    auto dim = h.rows();
    auto evec = std::vector<value_type>{};
    int verbosity = m_logger->max_trace_level == Logger::Info ? 3 : 0;
    itsolv::eigenproblem(evec, m_eigenvalues, h.data(), s.data(), dim, m_hermitian, m_svd_solver_threshold, verbosity);
    size_t n_solutions = 0;
    if (dim)
      n_solutions = evec.size() / dim;
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

  void solve_linear_equations(IXSpace<R, Q, P>& xspace) {
    m_logger->msg("SubspaceSolverLinEig::solve_linear_equations", Logger::Trace);
    auto h = xspace.data[EqnData::H];
    auto s = xspace.data[EqnData::S];
    auto rhs = xspace.data[EqnData::rhs];
    if (m_logger->data_dump) {
      m_logger->msg("S = " + as_string(s, 15), Logger::Info);
      m_logger->msg("H = " + as_string(h, 15), Logger::Info);
      m_logger->msg("rhs = " + as_string(rhs, 15), Logger::Info);
    }
    const auto dim = h.rows();
    const auto n_solutions = rhs.cols();
    auto solution = std::vector<value_type>{};
    m_eigenvalues.assign(n_solutions, 0);
    int verbosity = m_logger->max_trace_level == Logger::Info ? 3 : 0;
    itsolv::solve_LinearEquations(solution, m_eigenvalues, h.data(), s.data(), rhs.data(), dim, n_solutions,
                                  m_augmented_hessian, m_svd_solver_threshold, verbosity);
    m_solutions = Matrix<value_type>{std::move(solution), {n_solutions, dim}};
    m_errors.assign(size(), std::numeric_limits<value_type_abs>::max());
    if (m_logger->data_dump) {
      m_logger->msg("eigenvalues = ", begin(m_eigenvalues), end(m_eigenvalues), Logger::Debug, 10);
      m_logger->msg("solutions = " + as_string(m_solutions), Logger::Info);
    }
  }

public:
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

  // FIXME What difference does it make?
  //! Set Hermiticity of the subspace.
  void set_hermiticity(bool hermitian) { m_hermitian = hermitian; }
  bool get_hermiticity() { return m_hermitian; }
  //! Set value of augmented hessian parameter. If 0, than augmented Hessian is not used.
  void set_augmented_hessian(double parameter) { m_augmented_hessian = parameter; }
  double get_augmented_hessian() { return m_augmented_hessian; }

protected:
  Matrix<value_type> m_solutions;        //!< solution matrix with row vectors
  std::vector<value_type> m_eigenvalues; //!< eigenvalues
  std::vector<value_type_abs> m_errors;  //!< errors in subspace solutions
  std::shared_ptr<Logger> m_logger{};

public:
  value_type_abs m_svd_solver_threshold = 1.0e-14; //!< threshold to select null space during SVD in eigenproblem
protected:
  bool m_hermitian = false;       //!< flags the matrix as Hermitian
  double m_augmented_hessian = 0; //!< value of augmented hessian parameter. If 0, than augmented Hessian is not used
};

} // namespace molpro::linalg::itsolv::subspace

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SUBSPACESOLVERLINEIG_H
