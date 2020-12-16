#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SUBSPACESOLVEROPT_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SUBSPACESOLVEROPT_H
#include <molpro/linalg/itsolv/subspace/ISubspaceSolver.h>
#include <molpro/linalg/itsolv/subspace/IXSpace.h>
#include <molpro/linalg/itsolv/subspace/Matrix.h>

namespace molpro::linalg::itsolv::subspace {

/*!
 * @brief Solves subspace problem for linear eigenvalues and system of linear equations
 */
template <class RT, class QT, class PT>
class SubspaceSolverOpt : public ISubspaceSolver<RT, QT, PT> {
public:
  using value_type = typename ISubspaceSolver<RT, QT, PT>::value_type;
  using value_type_abs = typename ISubspaceSolver<RT, QT, PT>::value_type_abs;
  using R = typename ISubspaceSolver<RT, QT, PT>::R;
  using Q = typename ISubspaceSolver<RT, QT, PT>::Q;
  using P = typename ISubspaceSolver<RT, QT, PT>::P;

  explicit SubspaceSolverOpt(std::shared_ptr<Logger> logger) : m_logger(std::move(logger)) {}

  void solve(IXSpace<R, Q, P>& xspace, const size_t nroots_max) override {
    m_logger->msg("SubspaceSolverOpt::solve", Logger::Trace);
    assert(xspace.data.end() != xspace.data.find(EqnData::value));
    auto values = xspace.data[EqnData::value];
    assert(xspace.size() == values.size());

    if (true) {
      solve_steepest(xspace);
    } else {
    }
  }

protected:
  void solve_steepest(IXSpace<R, Q, P>& xspace) {
    m_logger->msg("SubspaceSolverOpt::solve_steepest", Logger::Trace);
    auto h = xspace.data[EqnData::H];
    auto s = xspace.data[EqnData::S];
    auto value = xspace.data[EqnData::value];
    if (m_logger->data_dump) {
      m_logger->msg("S = " + as_string(s), Logger::Info);
      m_logger->msg("H = " + as_string(h, 15), Logger::Info);
      m_logger->msg("value = " + as_string(value, 15), Logger::Info);
    }
    auto dim = h.rows();
    auto evec = std::vector<value_type>{};
    int verbosity = m_logger->max_trace_level == Logger::Info ? 3 : 0;
    m_solutions.resize({1, dim});
    m_solutions.slice().fill(0);
    m_solutions(0,dim-1)=1;
    m_errors.assign(1,h(0,0)); // FIXME
    if (m_logger->data_dump) {
      m_logger->msg("solution = " + as_string(m_solutions), Logger::Info);
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
  bool m_hermitian = false;        //!< flags the matrix as Hermitian
  double m_augmented_hessian = 0; //!< value of augmented hessian parameter. If 0, than augmented Hessian is not used
};

} // namespace molpro::linalg::itsolv::subspace

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SUBSPACESOLVEROPT_H
