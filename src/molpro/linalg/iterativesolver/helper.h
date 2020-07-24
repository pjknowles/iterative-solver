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
};

} // namespace iterativesolver
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_H_
