#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_H_
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_H_
#include <vector>

namespace molpro {
namespace linalg {
namespace iterativesolver {
template <typename scalar_type>
class helper {
public:
  static int propose_singularity_deletion(size_t n, size_t ndim, const scalar_type* m, const std::vector<int>& candidates,
                                   double threshold);
};
} // namespace iterativesolver
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_H_
