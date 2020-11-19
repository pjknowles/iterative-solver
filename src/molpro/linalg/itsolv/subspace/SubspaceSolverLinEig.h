#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SUBSPACESOLVERLINEIG_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SUBSPACESOLVERLINEIG_H
#include <molpro/linalg/itsolv/subspace/Matrix.h>
#include <molpro/linalg/itsolv/subspace/SubspaceSolverI.h>
#include <molpro/linalg/itsolv/subspace/XSpaceI.h>
#include <molpro/linalg/itsolv/subspace/check_conditioning.h>

namespace molpro::linalg::itsolv::subspace {
namespace detail {
// FIXME This can probably be done much better with eigen
/*!
 * @brief Transform square matrix to a new basis
 *
 * M' = L * (M * L.T)
 *
 * @param mat square matrix to transform
 * @param lin_trans linear transformation to a new basis, row vectors
 * @returns transformed matrix
 */
template <typename T>
auto transform(const Matrix<T>& mat, const Matrix<T>& lin_trans) {
  const size_t n = lin_trans.rows();
  const size_t m = lin_trans.cols();
  assert(mat.rows() == m && mat.cols() == m);
  auto mat1 = Matrix<T>({m, n});
  for (size_t i = 0; i < m; ++i)
    for (size_t j = 0; j < n; ++j)
      for (size_t k = 0; k < m; ++k)
        mat1(i, j) += mat(i, k) * lin_trans(j, k);
  auto result = Matrix<T>({n, n});
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j)
      for (size_t k = 0; k < m; ++k)
        result(i, j) += lin_trans(i, k) * mat1(k, j);
  return result;
}

/*!
 * @brief Transform solutions back to the original basis
 * Let solutions be {s}, transformed basis {t} and original basis {x}
 *
 * t_i = \sum_j L_{i,j} x_j
 * s_i = \sum_j C_{i,j} t_j
 *     = \sum_j C_{i,j} \sum_k L_{j,k} x_k
 *     = \sum_j C_{i,j} \sum_k L_{j,k} x_k
 *     = \sum_k K_{i,k} x_k
 *
 * K_{i,k} = \sum_j C_{i,j} L_{j,k}
 *
 * K = C L
 *
 * @param mat
 * @param lin_trans
 * @return
 */
template <typename T>
auto transform_solutions(const Matrix<T>& solutions, const Matrix<T>& lin_trans) {
  const size_t nsol = solutions.rows();
  const size_t n = lin_trans.rows();
  const size_t m = lin_trans.cols();
  assert(solutions.cols() == n);
  auto result = Matrix<T>({nsol, m});
  for (size_t i = 0; i < nsol; ++i)
    for (size_t j = 0; j < m; ++j)
      for (size_t k = 0; k < n; ++k)
        result(i, j) += solutions(i, k) * lin_trans(k, j);
  return result;
}

template <typename value_type, typename value_type_abs>
void normalise_transformation(Matrix<value_type>& m_lin_trans, const std::vector<value_type_abs>& norm,
                              value_type_abs threshold) {
  for (size_t i = 0; i < norm.size(); ++i) {
    if (norm[i] > threshold) {
      m_lin_trans.row(i).scal(1. / norm[i]);
    }
  }
}

} // namespace detail

template <class RT, class QT, class PT>
class SubspaceSolverLinEig : public SubspaceSolverI<RT, QT, PT> {
public:
  using value_type = typename SubspaceSolverI<RT, QT, PT>::value_type;
  using value_type_abs = typename SubspaceSolverI<RT, QT, PT>::value_type_abs;
  using R = typename SubspaceSolverI<RT, QT, PT>::R;
  using Q = typename SubspaceSolverI<RT, QT, PT>::Q;
  using P = typename SubspaceSolverI<RT, QT, PT>::P;

  explicit SubspaceSolverLinEig(std::shared_ptr<Logger> logger) : m_logger(std::move(logger)) {}

  void check_conditioning(XSpaceI<R, Q, P>& xspace) {
    auto nX_on_entry = xspace.dimensions().nX;
    auto nX = xspace.dimensions().nX;
    m_logger->msg("SubspaceSolverLinEig::check_conditioning size of x space = " + std::to_string(nX), Logger::Trace);
    if (m_logger->data_dump) {
      m_logger->msg("on entry", Logger::Info);
      m_logger->msg("Sxx = " + as_string(xspace.data[EqnData::S]), Logger::Info);
      m_logger->msg("Hxx = " + as_string(xspace.data[EqnData::H]), Logger::Info);
    }
    auto norm = xspace::check_conditioning_gram_schmidt(xspace, m_lin_trans, m_norm_stability_threshold, *m_logger);
    nX = xspace.dimensions().nX;
    m_logger->msg("size of x space after conditioning = " + std::to_string(nX), Logger::Debug);
    if (m_logger->data_dump) {
      auto n = Matrix<value_type_abs>{std::move(norm), {1, nX}};
      m_logger->msg("norm = " + as_string(n), Logger::Info);
    }
    if (m_logger->data_dump && nX != nX_on_entry) {
      m_logger->msg("Sxx = " + as_string(xspace.data[EqnData::S]), Logger::Info);
      m_logger->msg("Hxx = " + as_string(xspace.data[EqnData::H]), Logger::Info);
    }
  }

  void solve(XSpaceI<R, Q, P>& xspace, const size_t nroots_max) override {
    m_logger->msg("SubspaceSolverLinEig::solve", Logger::Trace);
    check_conditioning(xspace);
    auto h = detail::transform(xspace.data[EqnData::H], m_lin_trans);
    auto s = detail::transform(xspace.data[EqnData::S], m_lin_trans);
    if (m_logger->data_dump) {
      m_logger->msg("L = " + as_string(m_lin_trans), Logger::Info);
      m_logger->msg("S = " + as_string(s), Logger::Info);
      m_logger->msg("H = " + as_string(h), Logger::Info);
    }
    std::cout << "h "<<as_string(h)<<std::endl;
    std::cout << "s "<<as_string(s)<<std::endl;
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
    m_solutions = detail::transform_solutions(m_solutions, m_lin_trans);
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
  Matrix<value_type> m_lin_trans;        //!< linear transformation to a well conditioned subspace
  std::shared_ptr<Logger> m_logger;

public:
  value_type_abs m_norm_stability_threshold = 1.0e-4; //!< norm threshold for Gram Schmidt orthogonalisation
  value_type_abs m_svd_solver_threshold = 1.0e-14;    //!< threshold to select null space during SVD in eigenproblem
  bool m_hermitian = true;                            //!< flags the matrix as Hermitian
};

} // namespace molpro::linalg::itsolv::subspace

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SUBSPACESOLVERLINEIG_H
