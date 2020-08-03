#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_H_
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_H_
#include <complex>
#include <molpro/iostream.h>
#include <vector>

namespace molpro {
namespace linalg {
namespace iterativesolver {
template <typename T>
struct is_complex : std::false_type {};

template <typename T>
struct is_complex<std::complex<T>> : std::true_type {};

template <typename value_type>
int propose_singularity_deletion(size_t n, size_t ndim, const value_type* m, const std::vector<size_t>& candidates,
                                 double threshold);
template <typename value_type>
void printMatrix(const std::vector<value_type>, size_t rows, size_t cols, std::string title = "",
                 std::ostream& s = molpro::cout);

template <typename value_type, typename std::enable_if_t<is_complex<value_type>{}, int> = 0>
void eigenproblem(std::vector<value_type>& eigenvectors, std::vector<value_type>& eigenvalues,
                  const std::vector<value_type>& matrix, const std::vector<value_type>& metric, const size_t dimension,
                  bool hermitian, double svdThreshold, int verbosity = 0);

template <typename value_type, typename std::enable_if_t<!is_complex<value_type>{}, nullptr_t> = nullptr>
void eigenproblem(std::vector<value_type>& eigenvectors, std::vector<value_type>& eigenvalues,
                  const std::vector<value_type>& matrix, const std::vector<value_type>& metric, const size_t dimension,
                  bool hermitian, double svdThreshold, int verbosity = 0);

template <typename value_type, typename std::enable_if_t<is_complex<value_type>{}, int> = 0>
void solve_LinearEquations(std::vector<value_type>& solution, std::vector<value_type>& eigenvalues,
                           const std::vector<value_type>& matrix, const std::vector<value_type>& metric,
                           const std::vector<value_type> rhs, const size_t dimension, size_t nroot,
                           double augmented_hessian, double svdThreshold, int verbosity);

template <typename value_type, typename std::enable_if_t<!is_complex<value_type>{}, nullptr_t> = nullptr>
void solve_LinearEquations(std::vector<value_type>& solution, std::vector<value_type>& eigenvalues,
                           const std::vector<value_type>& matrix, const std::vector<value_type>& metric,
                           const std::vector<value_type> rhs, const size_t dimension, size_t nroot,
                           double augmented_hessian, double svdThreshold, int verbosity);

template <typename value_type>
void solve_DIIS(std::vector<value_type>& solution, const std::vector<value_type>& matrix, const size_t dimension,
                double svdThreshold, int verbosity = 0);

} // namespace iterativesolver
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_H_
