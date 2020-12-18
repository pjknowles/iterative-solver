#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SUBSPACESOLVEROPTBFGS_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SUBSPACESOLVEROPTBFGS_H
#include <molpro/linalg/itsolv/subspace/ISubspaceSolver.h>
#include <molpro/linalg/itsolv/subspace/IXSpace.h>
#include <molpro/linalg/itsolv/subspace/Matrix.h>

namespace molpro::linalg::itsolv::subspace {

/*!
 * @brief Solves subspace problem for minimisation using the Steepest Descent algorithm
 */
template <class RT, class QT, class PT>
class SubspaceSolverOptBFGS : public ISubspaceSolver<RT, QT, PT> {
public:
  using value_type = typename ISubspaceSolver<RT, QT, PT>::value_type;
  using value_type_abs = typename ISubspaceSolver<RT, QT, PT>::value_type_abs;
  using R = typename ISubspaceSolver<RT, QT, PT>::R;
  using Q = typename ISubspaceSolver<RT, QT, PT>::Q;
  using P = typename ISubspaceSolver<RT, QT, PT>::P;

  explicit SubspaceSolverOptBFGS(std::shared_ptr<Logger> logger) : m_logger(std::move(logger)) {}

  void solve(IXSpace<R, Q, P>& xspace, const size_t nroots_max) override {
    m_logger->msg("SubspaceSolverOptBFGS::solve", Logger::Trace);
    assert(xspace.data.end() != xspace.data.find(EqnData::value));
    auto values = xspace.data[EqnData::value];
    assert(xspace.size() == values.size());

    auto kH = xspace.data[EqnData::H];
    auto kS = xspace.data[EqnData::S];
    auto kValue = xspace.data[EqnData::value];
    if (m_logger->data_dump) {
      m_logger->msg("S = " + as_string(kS), Logger::Info);
      m_logger->msg("H = " + as_string(kH, 15), Logger::Info);
      m_logger->msg("value = " + as_string(kValue, 15), Logger::Info);
    }
    auto kDim = kH.rows();
    int kVerbosity = m_logger->max_trace_level == Logger::Info ? 3 : 0;
    m_solutions.resize({1, kDim});
    m_solutions.slice().fill(0);
    m_solutions(0, 0) = 1;
    m_errors.assign(1, kH(0, 0));
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
  const std::vector<value_type>& eigenvalues() const override {
    throw std::logic_error("eigenvalues() not available in non-linear method");
  }
  const std::vector<value_type_abs>& errors() const override { return m_errors; }

  //! Number of solutions
  size_t size() const override { return m_solutions.rows(); }

protected:
  Matrix<value_type> m_solutions;       //!< solution matrix with row vectors
  std::vector<value_type_abs> m_errors; //!< errors in subspace solutions
  std::shared_ptr<Logger> m_logger{};

public:
  value_type_abs m_svd_solver_threshold = 1.0e-14; //!< threshold to select null space during SVD in eigenproblem
};

} // namespace molpro::linalg::itsolv::subspace

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SUBSPACESOLVEROPTBFGS_H
