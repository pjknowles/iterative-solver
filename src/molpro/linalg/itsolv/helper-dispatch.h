#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_DISPATCH_H_
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_DISPATCH_H_

/*!
 * @file
 * @brief Precision-agnostic dispatch of the dense linear algebra kernels used by the subspace solvers.
 *
 * Every kernel in this file is templated on the scalar type and dispatches on it:
 *  - @c float, @c double, @c std::complex<float> and @c std::complex<double> are forwarded to LAPACK
 *    through the LAPACKE C interface, whenever one was found at build time (@c HAVE_LAPACKE);
 *  - any other scalar type is handled by the equivalent templated Eigen decomposition, which works for
 *    every scalar @c Eigen::NumTraits has been specialised for: @c long double out of the box, and
 *    extended- and arbitrary-precision types such as @c boost::multiprecision or @c mpfr::mpreal
 *    through the Eigen support those libraries provide.
 *
 * The two branches obey the same interface contract, so the caller never needs to know which one ran.
 */

#include <Eigen/Dense>

#include <molpro/lapacke.h>
#include <molpro/linalg/itsolv/helper.h>

#include <algorithm>
#include <complex>
#include <cstddef>
#include <list>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#ifdef MOLPRO
extern "C" int dsyev_c(char, char, int, double*, int, double*);
#endif

namespace molpro::linalg::itsolv {

/*!
 * @brief Whether the dense kernels of this library dispatch @p T to LAPACK rather than to Eigen.
 *
 * True for the four scalar types LAPACK itself provides kernels for, and only when an interface to
 * LAPACK is actually available in this build. Everything else -- in particular every extended- and
 * arbitrary-precision type -- goes through Eigen.
 */
template <typename T>
struct has_lapack_kernel : std::false_type {};
#if defined(HAVE_LAPACKE) || defined(MOLPRO)
template <>
struct has_lapack_kernel<double> : std::true_type {};
#endif
#ifdef HAVE_LAPACKE
template <>
struct has_lapack_kernel<float> : std::true_type {};
template <>
struct has_lapack_kernel<std::complex<float>> : std::true_type {};
template <>
struct has_lapack_kernel<std::complex<double>> : std::true_type {};
#endif
template <typename T>
inline constexpr bool has_lapack_kernel_v = has_lapack_kernel<T>::value;

#ifdef HAVE_LAPACKE
/*!
 * @brief Thin overload sets that hide the type-letter in the LAPACKE function names.
 *
 * std::complex<T> is layout-compatible with the corresponding LAPACKE complex type, whichever of the
 * several possible spellings (C99 _Complex, MKL_Complex8/16, a plain struct) the installed header uses.
 */
namespace lapack {

//! Eigen-decomposition of a real symmetric / complex hermitian matrix: divide and conquer (?syevd, ?heevd),
//! several times faster than the QR iteration of ?syev/?heev once the eigenvectors are wanted
inline lapack_int heev(int matrix_layout, char jobz, char uplo, lapack_int n, float* a, lapack_int lda, float* w) {
  return LAPACKE_ssyevd(matrix_layout, jobz, uplo, n, a, lda, w);
}
inline lapack_int heev(int matrix_layout, char jobz, char uplo, lapack_int n, double* a, lapack_int lda, double* w) {
  return LAPACKE_dsyevd(matrix_layout, jobz, uplo, n, a, lda, w);
}
inline lapack_int heev(int matrix_layout, char jobz, char uplo, lapack_int n, std::complex<float>* a, lapack_int lda,
                       float* w) {
  return LAPACKE_cheevd(matrix_layout, jobz, uplo, n, reinterpret_cast<lapack_complex_float*>(a), lda, w);
}
inline lapack_int heev(int matrix_layout, char jobz, char uplo, lapack_int n, std::complex<double>* a, lapack_int lda,
                       double* w) {
  return LAPACKE_zheevd(matrix_layout, jobz, uplo, n, reinterpret_cast<lapack_complex_double*>(a), lda, w);
}

//! Singular value decomposition, divide and conquer (?gesdd)
inline lapack_int gesdd(int matrix_layout, char jobz, lapack_int m, lapack_int n, float* a, lapack_int lda, float* s,
                        float* u, lapack_int ldu, float* vt, lapack_int ldvt) {
  return LAPACKE_sgesdd(matrix_layout, jobz, m, n, a, lda, s, u, ldu, vt, ldvt);
}
inline lapack_int gesdd(int matrix_layout, char jobz, lapack_int m, lapack_int n, double* a, lapack_int lda, double* s,
                        double* u, lapack_int ldu, double* vt, lapack_int ldvt) {
  return LAPACKE_dgesdd(matrix_layout, jobz, m, n, a, lda, s, u, ldu, vt, ldvt);
}
inline lapack_int gesdd(int matrix_layout, char jobz, lapack_int m, lapack_int n, std::complex<float>* a,
                        lapack_int lda, float* s, std::complex<float>* u, lapack_int ldu, std::complex<float>* vt,
                        lapack_int ldvt) {
  return LAPACKE_cgesdd(matrix_layout, jobz, m, n, reinterpret_cast<lapack_complex_float*>(a), lda, s,
                        reinterpret_cast<lapack_complex_float*>(u), ldu, reinterpret_cast<lapack_complex_float*>(vt),
                        ldvt);
}
inline lapack_int gesdd(int matrix_layout, char jobz, lapack_int m, lapack_int n, std::complex<double>* a,
                        lapack_int lda, double* s, std::complex<double>* u, lapack_int ldu, std::complex<double>* vt,
                        lapack_int ldvt) {
  return LAPACKE_zgesdd(matrix_layout, jobz, m, n, reinterpret_cast<lapack_complex_double*>(a), lda, s,
                        reinterpret_cast<lapack_complex_double*>(u), ldu, reinterpret_cast<lapack_complex_double*>(vt),
                        ldvt);
}

//! Singular value decomposition, QR iteration (?gesvd)
inline lapack_int gesvd(int matrix_layout, char jobu, char jobvt, lapack_int m, lapack_int n, float* a, lapack_int lda,
                        float* s, float* u, lapack_int ldu, float* vt, lapack_int ldvt, float* superb) {
  return LAPACKE_sgesvd(matrix_layout, jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, superb);
}
inline lapack_int gesvd(int matrix_layout, char jobu, char jobvt, lapack_int m, lapack_int n, double* a, lapack_int lda,
                        double* s, double* u, lapack_int ldu, double* vt, lapack_int ldvt, double* superb) {
  return LAPACKE_dgesvd(matrix_layout, jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, superb);
}
inline lapack_int gesvd(int matrix_layout, char jobu, char jobvt, lapack_int m, lapack_int n, std::complex<float>* a,
                        lapack_int lda, float* s, std::complex<float>* u, lapack_int ldu, std::complex<float>* vt,
                        lapack_int ldvt, float* superb) {
  return LAPACKE_cgesvd(matrix_layout, jobu, jobvt, m, n, reinterpret_cast<lapack_complex_float*>(a), lda, s,
                        reinterpret_cast<lapack_complex_float*>(u), ldu, reinterpret_cast<lapack_complex_float*>(vt),
                        ldvt, superb);
}
inline lapack_int gesvd(int matrix_layout, char jobu, char jobvt, lapack_int m, lapack_int n, std::complex<double>* a,
                        lapack_int lda, double* s, std::complex<double>* u, lapack_int ldu, std::complex<double>* vt,
                        lapack_int ldvt, double* superb) {
  return LAPACKE_zgesvd(matrix_layout, jobu, jobvt, m, n, reinterpret_cast<lapack_complex_double*>(a), lda, s,
                        reinterpret_cast<lapack_complex_double*>(u), ldu, reinterpret_cast<lapack_complex_double*>(vt),
                        ldvt, superb);
}

} // namespace lapack
#endif // HAVE_LAPACKE

namespace detail {

/*!
 * @brief Eigen implementation of the hermitian eigenproblem, used for every scalar type LAPACK does
 * not cover.
 *
 * @param[in,out] a on entry the matrix in column-major order, on exit the eigenvectors as its columns
 * @param[out] w the eigenvalues in ascending order
 * @param[in] dimension length of one axis of the matrix
 * \returns 0 on success, 1 if the decomposition did not converge
 */
template <typename value_type>
int eigensolver_hermitian_kernel(std::false_type /*use_lapack*/, std::span<value_type> a,
                                 std::span<real_type_t<value_type>> w, size_t dimension) {
  using matrix_type = Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
  using real_vector_type = Eigen::Vector<real_type_t<value_type>, Eigen::Dynamic>;
  Eigen::Map<matrix_type> A(a.data(), dimension, dimension);
  // Eigen, like ?syev with uplo='L', references only the lower triangle and returns ascending eigenvalues
  Eigen::SelfAdjointEigenSolver<matrix_type> solver(A, Eigen::ComputeEigenvectors);
  if (solver.info() != Eigen::Success)
    return 1;
  Eigen::Map<real_vector_type>(w.data(), dimension) = solver.eigenvalues();
  A = solver.eigenvectors();
  return 0;
}

#if defined(HAVE_LAPACKE) || defined(MOLPRO)
//! LAPACK implementation of the hermitian eigenproblem, @see the Eigen overload for the contract
template <typename value_type>
int eigensolver_hermitian_kernel(std::true_type /*use_lapack*/, std::span<value_type> a,
                                 std::span<real_type_t<value_type>> w, size_t dimension) {
  constexpr char compute_eigenvalues_eigenvectors = 'V';
  constexpr char store_lower_triangle = 'L';
#ifdef MOLPRO
  if constexpr (std::is_same_v<value_type, double>) {
    return dsyev_c(compute_eigenvalues_eigenvectors, store_lower_triangle, int(dimension), a.data(), int(dimension),
                   w.data());
  }
#endif
#ifdef HAVE_LAPACKE
  return int(lapack::heev(LAPACK_COL_MAJOR, compute_eigenvalues_eigenvectors, store_lower_triangle,
                          lapack_int(dimension), a.data(), lapack_int(dimension), w.data()));
#else
  throw std::logic_error("no LAPACK kernel available for this scalar type");
#endif
}
#endif

} // namespace detail

/*!
 * @brief Eigen-decomposition of a real symmetric / complex hermitian matrix.
 *
 * Dispatches to ?syevd/?heevd for the scalar types LAPACK covers and to Eigen::SelfAdjointEigenSolver
 * for every other precision; both branches produce the same layout and ordering.
 *
 * @param[in] matrix the input matrix (will not be altered). Must be square, dimension*dimension elements long.
 *                   Only the lower triangle in column-major order is referenced.
 * @param[out] eigenvectors the eigenvectors as the columns of a column-major matrix, i.e. eigenvector
 *                   @c i occupies elements @c [i*dimension, (i+1)*dimension). Same size as matrix.
 * @param[out] eigenvalues the eigenvalues in ascending order. Must be @p dimension elements long.
 * @param[in] dimension length of one axis of the matrix.
 * \returns status. If 0, successful exit. If -i, the ith argument had an illegal value. If i, the
 *          algorithm failed to converge.
 */
template <typename value_type>
int eigensolver_hermitian(std::span<const value_type> matrix, std::span<value_type> eigenvectors,
                          std::span<real_type_t<value_type>> eigenvalues, const size_t dimension) {
  // validate input
  if (eigenvectors.size() != matrix.size()) {
    throw std::runtime_error("Matrix of eigenvectors and input matrix are not the same size! (" +
                             std::to_string(eigenvectors.size()) + " vs. " + std::to_string(matrix.size()) + ")");
  }

  if (eigenvectors.size() != dimension * dimension || eigenvalues.size() != dimension) {
    throw std::runtime_error("Size of eigenvectors/eigenvalues do not match dimension!");
  }

  // copy input matrix, since the decomposition overwrites it in place
  std::copy(matrix.begin(), matrix.end(), eigenvectors.begin());
  if (dimension == 0)
    return 0;

  return detail::eigensolver_hermitian_kernel<value_type>(has_lapack_kernel<value_type>{}, eigenvectors, eigenvalues,
                                                          dimension);
}

/*!
 * @brief Eigen-decomposition of a real symmetric / complex hermitian matrix, @see eigensolver_hermitian.
 *
 * @param[in] dimension length of one axis of the matrix.
 * @param[in] matrix the input matrix (will not be altered). Must be square, dimension*dimension elements long.
 * \returns a std::list of instances of SVD, a struct containing one eigenvalue and one eigenvector,
 *          in descending order of the eigenvalue. For a hermitian positive semi-definite matrix these
 *          are the singular values and right singular vectors.
 */
template <typename value_type>
std::list<SVD<value_type>> eigensolver_hermitian(size_t dimension, std::span<const value_type> matrix) {
  using real_t = real_type_t<value_type>;
  std::vector<value_type> eigvecs(dimension * dimension);
  std::vector<real_t> eigvals(dimension);

  const int success = eigensolver_hermitian<value_type>(matrix, eigvecs, eigvals, dimension);
  if (success < 0) {
    throw std::invalid_argument("Invalid argument of eigensolver_hermitian: " + std::to_string(-success));
  }
  if (success > 0) {
    throw std::runtime_error("Hermitian eigensolver failed to converge. "
                             " elements of an intermediate tridiagonal form did not converge to zero.");
  }

  auto eigensystem = std::list<SVD<value_type>>{};

  // populate eigensystem
  for (int i = int(dimension) - 1; i >= 0;
       i--) { // note: flipping this axis gives parity with results of eigen::jacobiSVD
    auto temp_eigenproblem = SVD<value_type>{};
    temp_eigenproblem.value = eigvals[i];
    for (size_t j = 0; j < dimension; j++) {
      temp_eigenproblem.v.emplace_back(eigvecs[j + (dimension * i)]);
    }
    eigensystem.emplace_back(temp_eigenproblem);
  }

  return eigensystem;
}

} // namespace molpro::linalg::itsolv

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_DISPATCH_H_
