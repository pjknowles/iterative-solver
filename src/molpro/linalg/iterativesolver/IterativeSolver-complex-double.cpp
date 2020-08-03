#include <complex>
#include <molpro/linalg/iterativesolver/helper-implementation.h>
namespace {
using value_type = std::complex<double>;
}
namespace molpro {
namespace linalg {
namespace iterativesolver {

template int propose_singularity_deletion<value_type>(size_t n, size_t ndim, const value_type* m,
                                                      const std::vector<size_t>& candidates, double threshold);

template void printMatrix<value_type>(const std::vector<value_type>, size_t rows, size_t cols, std::string title = "",
                                      std::ostream& s = molpro::cout);

template void eigenproblem<value_type>(std::vector<value_type>& eigenvectors, std::vector<value_type>& eigenvalues,
                                       const std::vector<value_type>& matrix, const std::vector<value_type>& metric,
                                       const size_t dimension, bool hermitian, double svdThreshold, int verbosity = 0);

template void solve_LinearEquations<value_type>(std::vector<value_type>& solution, std::vector<value_type>& eigenvalues,
                                                const std::vector<value_type>& matrix,
                                                const std::vector<value_type>& metric,
                                                const std::vector<value_type> rhs, const size_t dimension, size_t nroot,
                                                double augmented_hessian, double svdThreshold, int verbosity);

template void solve_DIIS<value_type>(std::vector<value_type>& solution, const std::vector<value_type>& matrix,
                                     const size_t dimension, double svdThreshold, int verbosity = 0);
} // namespace iterativesolver
} // namespace linalg
} // namespace molpro
