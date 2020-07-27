#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_H_
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_H_
#include <molpro/iostream.h>
#include <vector>

namespace molpro {
namespace linalg {
namespace iterativesolver {
template <typename value_type>
class helper {
public:
  static int propose_singularity_deletion(size_t n, size_t ndim, const value_type* m,
                                          const std::vector<int>& candidates, double threshold);
  static void printMatrix(const std::vector<value_type>, size_t rows, size_t cols, std::string title = "",
                          std::ostream& s = molpro::cout);

  static void eigenproblem(std::vector<value_type>& eigenvectors, std::vector<value_type>& eigenvalues,
                           const std::vector<value_type>& matrix, const std::vector<value_type>& metric,
                           const size_t dimension, bool hermitian, double svdThreshold, int verbosity = 0);

  static void solve_DIIS(std::vector<value_type>& solution, const std::vector<value_type>& matrix,
                         const size_t nQ, double svdThreshold, int verbosity=0);
};

} // namespace iterativesolver
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_H_
