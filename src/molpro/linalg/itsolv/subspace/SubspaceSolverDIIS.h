#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SUBSPACESOLVERDIIS_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SUBSPACESOLVERDIIS_H
#include <algorithm>
#include <molpro/linalg/itsolv/subspace/ISubspaceSolver.h>
#include <molpro/linalg/itsolv/subspace/IXSpace.h>
#include <molpro/linalg/itsolv/subspace/Matrix.h>

namespace molpro::linalg::itsolv::subspace {

/*!
 * @brief Solves non-linear equations using DIIS
 */
template <class RT, class QT, class PT>
class SubspaceSolverDIIS : public ISubspaceSolver<RT, QT, PT> {
  const bool& m_converged;

public:
  using value_type = typename ISubspaceSolver<RT, QT, PT>::value_type;
  using value_type_abs = typename ISubspaceSolver<RT, QT, PT>::value_type_abs;
  using R = typename ISubspaceSolver<RT, QT, PT>::R;
  using Q = typename ISubspaceSolver<RT, QT, PT>::Q;
  using P = typename ISubspaceSolver<RT, QT, PT>::P;

  explicit SubspaceSolverDIIS(std::shared_ptr<Logger> logger, const bool& converged)
      : m_converged(converged), m_logger(std::move(logger)) {}

  void solve(IXSpace<R, Q, P>& xspace, const size_t nroots_max) override {
    m_logger->msg("SubspaceSolverDIIS::solve", Logger::Trace);

    auto kH = xspace.data[EqnData::H];
    auto kS = xspace.data[EqnData::S];
    if (m_logger->data_dump) {
      m_logger->msg("S = " + as_string(kS), Logger::Info);
      m_logger->msg("H = " + as_string(kH, 15), Logger::Info);
    }
    auto kDim = kH.rows();
    m_solutions.resize({1, kDim});
    if (m_converged) {
      m_solutions.fill(0);
      m_solutions(0, 0) = 1;
      return;
    }
    //    int kVerbosity = m_logger->max_trace_level == Logger::Info ? 3 : 0;
    //    auto dH = kH;
    //    auto dDim = kDim - 1;
    //    dH.resize({dDim, dDim});
    //    for (size_t i = 0; i < dDim; i++)
    //      for (size_t j = 0; j < dDim; j++)
    //        dH(i, j) = kH(i, j) - kH(i + 1, j) - kH(i, j + 1) + kH(i + 1, j + 1);
    //    m_logger->msg("dH = " + as_string(kH, 15), Logger::Info);
    std::vector<value_type> solution(kDim);
    std::vector<value_type> matrix;
    matrix.reserve(kDim * kDim);
    //    std::copy(std::begin(kH),std::end(kH),matrix.begin());
    for (size_t i = 0; i < kDim; ++i)
      for (size_t j = 0; j < kDim; ++j)
        matrix.push_back(kH(j, i));
    solve_DIIS(solution, matrix, kDim, 1e-10, 1);
    //    std::copy(solution.begin(),solution.end(),m_solutions.begin());
    for (size_t i = 0; i < kDim; ++i)
      m_solutions(0, i) = solution[i];
    m_errors.assign(1, kH(0, 0)); // TODO fix
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

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SUBSPACESOLVERDIIS_H
