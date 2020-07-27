#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_H_
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_H_
#include <vector>
#include <molpro/iostream.h>

namespace molpro {
namespace linalg {
namespace iterativesolver {
template <typename scalar_type>
class helper {
public:
  static int propose_singularity_deletion(size_t n, size_t ndim, const scalar_type* m, const std::vector<int>& candidates,
                                   double threshold);
  static void printMatrix(const std::vector<scalar_type>, size_t rows, size_t cols, std::string title="", std::ostream& s=molpro::cout);

  static void eigenproblem(std::vector<scalar_type>& eigenvectors, std::vector<scalar_type>& eigenvalues,const std::vector<scalar_type> matrix, const std::vector<scalar_type> metric, const size_t dimension,  bool hermitian, double svdThreshold, int verbosity=0);
};

} // namespace iterativesolver
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_H_
