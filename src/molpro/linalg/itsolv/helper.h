#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_H_
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_H_
#include <complex>
#include <cstddef>
#include <list>
#include <molpro/iostream.h>
#include <molpro/linalg/array/Span.h>
#include <vector>

namespace molpro {
namespace linalg {
namespace itsolv {
template <typename T>
struct is_complex : std::false_type {};

template <typename T>
struct is_complex<std::complex<T>> : std::true_type {};

template <typename value_type>
int propose_singularity_deletion(size_t n, size_t ndim, const value_type* m, const std::vector<size_t>& candidates,
                                 double threshold);

//! Stores a singular value and corresponding left and right singular vectors
template <typename T>
struct SVD {
  using value_type = T;
  value_type value;
  std::vector<value_type> u; //!< left singular vector
  std::vector<value_type> v; //!< right singular vector
};

/*!
 * @brief Performs singular value decomposition and returns SVD objects for singular values less than threshold, sorted
 * in ascending order
 * @tparam value_type
 * @param ndim dimension size of a square matrix
 * @param m row-wise data buffer for a square matrix
 * @param threshold singular values less than threshold will be returned
 */
template <typename value_type, typename std::enable_if_t<!is_complex<value_type>{}, std::nullptr_t> = nullptr>
std::list<SVD<value_type>> svd_system(size_t ndim, const array::Span<value_type>& m, double threshold);
template <typename value_type, typename std::enable_if_t<is_complex<value_type>{}, int> = 0>
std::list<SVD<value_type>> svd_system(size_t ndim, const array::Span<value_type>& m, double threshold);

template <typename value_type>
void printMatrix(const std::vector<value_type>&, size_t rows, size_t cols, std::string title = "",
                 std::ostream& s = molpro::cout);

template <typename value_type, typename std::enable_if_t<is_complex<value_type>{}, int> = 0>
void eigenproblem(std::vector<value_type>& eigenvectors, std::vector<value_type>& eigenvalues,
                  const std::vector<value_type>& matrix, const std::vector<value_type>& metric, size_t dimension,
                  bool hermitian, double svdThreshold, int verbosity = 0);

template <typename value_type, typename std::enable_if_t<!is_complex<value_type>{}, std::nullptr_t> = nullptr>
void eigenproblem(std::vector<value_type>& eigenvectors, std::vector<value_type>& eigenvalues,
                  const std::vector<value_type>& matrix, const std::vector<value_type>& metric, size_t dimension,
                  bool hermitian, double svdThreshold, int verbosity = 0);

template <typename value_type, typename std::enable_if_t<is_complex<value_type>{}, int> = 0>
void solve_LinearEquations(std::vector<value_type>& solution, std::vector<value_type>& eigenvalues,
                           const std::vector<value_type>& matrix, const std::vector<value_type>& metric,
                           const std::vector<value_type>& rhs, size_t dimension, size_t nroot, double augmented_hessian,
                           double svdThreshold, int verbosity);

template <typename value_type, typename std::enable_if_t<!is_complex<value_type>{}, std::nullptr_t> = nullptr>
void solve_LinearEquations(std::vector<value_type>& solution, std::vector<value_type>& eigenvalues,
                           const std::vector<value_type>& matrix, const std::vector<value_type>& metric,
                           const std::vector<value_type>& rhs, size_t dimension, size_t nroot, double augmented_hessian,
                           double svdThreshold, int verbosity);

template <typename value_type>
void solve_DIIS(std::vector<value_type>& solution, const std::vector<value_type>& matrix, size_t dimension,
                double svdThreshold, int verbosity = 0);

/*
 * Explicit instantiation of double type
 */
extern template int propose_singularity_deletion<double>(size_t n, size_t ndim, const double* m,
                                                         const std::vector<size_t>& candidates, double threshold);

extern template void printMatrix<double>(const std::vector<double>&, size_t rows, size_t cols, std::string title,
                                         std::ostream& s);

extern template std::list<SVD<double>> svd_system(size_t ndim, const array::Span<double>& m, double threshold);

extern template void eigenproblem<double>(std::vector<double>& eigenvectors, std::vector<double>& eigenvalues,
                                          const std::vector<double>& matrix, const std::vector<double>& metric,
                                          const size_t dimension, bool hermitian, double svdThreshold,
                                          int verbosity = 0);

extern template void eigenproblem<double>(std::vector<double>& eigenvectors, std::vector<double>& eigenvalues,
                                          const std::vector<double>& matrix, const std::vector<double>& metric,
                                          const size_t dimension, bool hermitian, double svdThreshold,
                                          int verbosity = 0);

extern template void solve_LinearEquations<double>(std::vector<double>& solution, std::vector<double>& eigenvalues,
                                                   const std::vector<double>& matrix, const std::vector<double>& metric,
                                                   const std::vector<double>& rhs, size_t dimension, size_t nroot,
                                                   double augmented_hessian, double svdThreshold, int verbosity);

extern template void solve_LinearEquations<double>(std::vector<double>& solution, std::vector<double>& eigenvalues,
                                                   const std::vector<double>& matrix, const std::vector<double>& metric,
                                                   const std::vector<double>& rhs, size_t dimension, size_t nroot,
                                                   double augmented_hessian, double svdThreshold, int verbosity);

extern template void solve_DIIS<double>(std::vector<double>& solution, const std::vector<double>& matrix,
                                        const size_t dimension, double svdThreshold, int verbosity = 0);

/*
 * Explicit instantiation of std::complex<double> type
 */
extern template int propose_singularity_deletion<std::complex<double>>(size_t n, size_t ndim,
                                                                       const std::complex<double>* m,
                                                                       const std::vector<size_t>& candidates,
                                                                       double threshold);

extern template void printMatrix<std::complex<double>>(const std::vector<std::complex<double>>&, size_t rows,
                                                       size_t cols, std::string title, std::ostream& s);

extern template std::list<SVD<std::complex<double>>> svd_system(size_t ndim, const array::Span<std::complex<double>>& m,
                                                                double threshold);

extern template void eigenproblem<std::complex<double>>(std::vector<std::complex<double>>& eigenvectors,
                                                        std::vector<std::complex<double>>& eigenvalues,
                                                        const std::vector<std::complex<double>>& matrix,
                                                        const std::vector<std::complex<double>>& metric,
                                                        const size_t dimension, bool hermitian, double svdThreshold,
                                                        int verbosity = 0);

extern template void eigenproblem<std::complex<double>>(std::vector<std::complex<double>>& eigenvectors,
                                                        std::vector<std::complex<double>>& eigenvalues,
                                                        const std::vector<std::complex<double>>& matrix,
                                                        const std::vector<std::complex<double>>& metric,
                                                        const size_t dimension, bool hermitian, double svdThreshold,
                                                        int verbosity = 0);

extern template void solve_LinearEquations<std::complex<double>>(
    std::vector<std::complex<double>>& solution, std::vector<std::complex<double>>& eigenvalues,
    const std::vector<std::complex<double>>& matrix, const std::vector<std::complex<double>>& metric,
    const std::vector<std::complex<double>>& rhs, size_t dimension, size_t nroot, double augmented_hessian,
    double svdThreshold, int verbosity);

extern template void solve_LinearEquations<std::complex<double>>(
    std::vector<std::complex<double>>& solution, std::vector<std::complex<double>>& eigenvalues,
    const std::vector<std::complex<double>>& matrix, const std::vector<std::complex<double>>& metric,
    const std::vector<std::complex<double>>& rhs, size_t dimension, size_t nroot, double augmented_hessian,
    double svdThreshold, int verbosity);

extern template void solve_DIIS<std::complex<double>>(std::vector<std::complex<double>>& solution,
                                                      const std::vector<std::complex<double>>& matrix,
                                                      const size_t dimension, double svdThreshold, int verbosity = 0);
} // namespace itsolv
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_H_
