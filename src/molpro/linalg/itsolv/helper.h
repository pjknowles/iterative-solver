#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_H_
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_H_
#include <complex>
#include <cstddef>
#include <limits>
#include <list>
#include <molpro/iostream.h>
#include <molpro/linalg/array/Span.h>
#include <molpro/linalg/scalar_traits.h>
#include <span>
#include <type_traits>
#include <vector>

namespace molpro::linalg::itsolv {

// Scalar traits shared with molpro::linalg::array, re-exported here for backwards compatibility
using molpro::linalg::conjugate;
using molpro::linalg::is_complex;
using molpro::linalg::imaginary_part;
using molpro::linalg::precision_scaled;
using molpro::linalg::real_part;
using molpro::linalg::real_type;
using molpro::linalg::real_type_t;

//! Stores a singular value and corresponding left and right singular vectors
template <typename T>
struct SVD {
  using value_type = T;
  value_type value;
  std::vector<value_type> u; //!< left singular vector
  std::vector<value_type> v; //!< right singular vector
};

/*!
 * @brief Eigen-decomposition of a real symmetric matrix in double precision.
 *
 * Retained for backwards compatibility; equivalent to eigensolver_hermitian<double>(), which is the
 * precision-generic entry point. Note that despite the name the decomposition is only carried out by
 * LAPACK when this build found an interface to it.
 */
int eigensolver_lapacke_dsyev(std::span<const double> matrix, std::span<double> eigenvectors,
                              std::span<double> eigenvalues, const size_t dimension);

//! @copydoc eigensolver_lapacke_dsyev
std::list<SVD<double>> eigensolver_lapacke_dsyev(size_t dimension, std::span<const double> matrix);

template <typename value_type>
size_t get_rank(std::span<const value_type> eigenvalues, value_type threshold);

template <typename value_type>
size_t get_rank(std::list<SVD<value_type>> svd_system, real_type_t<value_type> threshold);

//! @copydoc get_rank
template <typename value_type>
size_t get_rank(const std::vector<value_type>& eigenvalues, value_type threshold) {
  return get_rank<value_type>(std::span<const value_type>{eigenvalues.data(), eigenvalues.size()}, threshold);
}

/*!
 * @brief Performs singular value decomposition and returns SVD objects for singular values less than threshold, sorted
 * in ascending order
 * @tparam value_type
 * @param nrows number of rows in the matrix
 * @param ncols number of columns in the matrix
 * @param m row-wise data buffer for a matrix
 * @param threshold singular values less than threshold will be returned
 */
template <typename value_type, typename std::enable_if_t<!is_complex<value_type>{}, std::nullptr_t> = nullptr>
std::list<SVD<value_type>> svd_system(size_t nrows, size_t ncols, const array::Span<value_type>& m,
                                      real_type_t<value_type> threshold, bool hermitian = false,
                                      bool reduce_to_rank = false);
template <typename value_type, typename std::enable_if_t<is_complex<value_type>{}, int> = 0>
std::list<SVD<value_type>> svd_system(size_t nrows, size_t ncols, const array::Span<value_type>& m,
                                      real_type_t<value_type> threshold, bool hermitian = false,
                                      bool reduce_to_rank = false);

template <typename value_type>
void printMatrix(const std::vector<value_type>&, size_t rows, size_t cols, std::string title = "",
                 std::ostream& s = molpro::cout);

template <typename value_type, typename std::enable_if_t<is_complex<value_type>{}, int> = 0>
void eigenproblem(std::vector<value_type>& eigenvectors, std::vector<value_type>& eigenvalues,
                  const std::vector<value_type>& matrix, const std::vector<value_type>& metric, size_t dimension,
                  bool hermitian, real_type_t<value_type> svdThreshold, int verbosity,
                  std::vector<std::pair<std::size_t, value_type>>* imag_eval_parts = nullptr,
                  size_t nroots = std::numeric_limits<size_t>::max());

/*!
 * @param imag_eval_parts if not null, receives a (root, imaginary part) record for every eigenvalue that is a member of a
 * complex conjugate pair. The members of a pair are adjacent roots and their records are consecutive, first member
 * first. Their eigenvalues are set to the common real part and their eigenvectors to the real and imaginary part of the
 * first member's complex eigenvector.
 */
template <typename value_type, typename std::enable_if_t<!is_complex<value_type>{}, std::nullptr_t> = nullptr>
void eigenproblem(std::vector<value_type>& eigenvectors, std::vector<value_type>& eigenvalues,
                  const std::vector<value_type>& matrix, const std::vector<value_type>& metric, size_t dimension,
                  bool hermitian, real_type_t<value_type> svdThreshold, int verbosity,
                  std::vector<std::pair<std::size_t, value_type>>* imag_eval_parts = nullptr,
                  size_t nroots = std::numeric_limits<size_t>::max());

template <typename value_type, typename std::enable_if_t<is_complex<value_type>{}, int> = 0>
void solve_LinearEquations(std::vector<value_type>& solution, std::vector<value_type>& eigenvalues,
                           const std::vector<value_type>& matrix, const std::vector<value_type>& metric,
                           const std::vector<value_type>& rhs, size_t dimension, size_t nroot,
                           real_type_t<value_type> augmented_hessian, real_type_t<value_type> svdThreshold,
                           int verbosity);

template <typename value_type, typename std::enable_if_t<!is_complex<value_type>{}, std::nullptr_t> = nullptr>
void solve_LinearEquations(std::vector<value_type>& solution, std::vector<value_type>& eigenvalues,
                           const std::vector<value_type>& matrix, const std::vector<value_type>& metric,
                           const std::vector<value_type>& rhs, size_t dimension, size_t nroot,
                           real_type_t<value_type> augmented_hessian, real_type_t<value_type> svdThreshold,
                           int verbosity);

template <typename value_type, typename std::enable_if_t<is_complex<value_type>{}, int> = 0>
void solve_DIIS(std::vector<value_type>& solution, const std::vector<value_type>& matrix, size_t dimension,
                real_type_t<value_type> svdThreshold, int verbosity = 0);
template <typename value_type, typename std::enable_if_t<!is_complex<value_type>{}, std::nullptr_t> = nullptr>
void solve_DIIS(std::vector<value_type>& solution, const std::vector<value_type>& matrix, size_t dimension,
                real_type_t<value_type> svdThreshold, int verbosity = 0);

/*
 * Explicit instantiation of double type
 */

extern template void printMatrix<double>(const std::vector<double>&, size_t rows, size_t cols, std::string title,
                                         std::ostream& s);

extern template size_t get_rank<double>(std::span<const double> eigenvalues, double threshold);

extern template std::list<SVD<double>> svd_system(size_t nrows, size_t ncols, const array::Span<double>& m,
                                                  double threshold, bool hermitian, bool reduce_to_rank);

extern template void eigenproblem<double>(std::vector<double>& eigenvectors, std::vector<double>& eigenvalues,
                                          const std::vector<double>& matrix, const std::vector<double>& metric,
                                          const size_t dimension, bool hermitian, double svdThreshold, int verbosity,
                                          std::vector<std::pair<std::size_t, double>> *imag_eval_parts, size_t nroots);

extern template void solve_LinearEquations<double>(std::vector<double>& solution, std::vector<double>& eigenvalues,
                                                   const std::vector<double>& matrix, const std::vector<double>& metric,
                                                   const std::vector<double>& rhs, size_t dimension, size_t nroot,
                                                   double augmented_hessian, double svdThreshold, int verbosity);

extern template void solve_DIIS<double>(std::vector<double>& solution, const std::vector<double>& matrix,
                                        const size_t dimension, double svdThreshold, int verbosity);

/*
 * Explicit instantiation of long double type
 *
 * The dense kernels fall back on Eigen for this precision, as they do for any other scalar type for
 * which LAPACK provides no kernel.
 */

extern template void printMatrix<long double>(const std::vector<long double>&, size_t rows, size_t cols,
                                              std::string title, std::ostream& s);

extern template size_t get_rank<long double>(std::span<const long double> eigenvalues, long double threshold);

extern template std::list<SVD<long double>> svd_system(size_t nrows, size_t ncols, const array::Span<long double>& m,
                                                       long double threshold, bool hermitian, bool reduce_to_rank);

extern template void eigenproblem<long double>(std::vector<long double>& eigenvectors,
                                               std::vector<long double>& eigenvalues,
                                               const std::vector<long double>& matrix,
                                               const std::vector<long double>& metric, const size_t dimension,
                                               bool hermitian, long double svdThreshold, int verbosity,
                                               std::vector<std::pair<std::size_t, long double>> *imag_eval_parts, size_t nroots);

extern template void solve_LinearEquations<long double>(
    std::vector<long double>& solution, std::vector<long double>& eigenvalues, const std::vector<long double>& matrix,
    const std::vector<long double>& metric, const std::vector<long double>& rhs, size_t dimension, size_t nroot,
    long double augmented_hessian, long double svdThreshold, int verbosity);

extern template void solve_DIIS<long double>(std::vector<long double>& solution, const std::vector<long double>& matrix,
                                             const size_t dimension, long double svdThreshold, int verbosity);

/*
 * Explicit instantiation of std::complex<double> type
 */
extern template void printMatrix<std::complex<double>>(const std::vector<std::complex<double>>&, size_t rows,
                                                       size_t cols, std::string title, std::ostream& s);

extern template std::list<SVD<std::complex<double>>> svd_system(size_t nrows, size_t ncols,
                                                                const array::Span<std::complex<double>>& m,
                                                                double threshold, bool hermitian, bool reduce_to_rank);

extern template void eigenproblem<std::complex<double>>(
    std::vector<std::complex<double>>& eigenvectors, std::vector<std::complex<double>>& eigenvalues,
    const std::vector<std::complex<double>>& matrix, const std::vector<std::complex<double>>& metric,
    const size_t dimension, bool hermitian, double svdThreshold, int verbosity,
    std::vector<std::pair<std::size_t, std::complex<double>>>* imag_eval_parts, size_t nroots);

extern template void solve_LinearEquations<std::complex<double>>(
    std::vector<std::complex<double>>& solution, std::vector<std::complex<double>>& eigenvalues,
    const std::vector<std::complex<double>>& matrix, const std::vector<std::complex<double>>& metric,
    const std::vector<std::complex<double>>& rhs, size_t dimension, size_t nroot, double augmented_hessian,
    double svdThreshold, int verbosity);

extern template void solve_DIIS<std::complex<double>>(std::vector<std::complex<double>>& solution,
                                                      const std::vector<std::complex<double>>& matrix,
                                                      const size_t dimension, double svdThreshold, int verbosity);
} // namespace molpro::linalg::itsolv
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_H_
