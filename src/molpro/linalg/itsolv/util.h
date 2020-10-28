#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_UTIL_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_UTIL_H
#include <molpro/linalg/itsolv/subspace/Matrix.h>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace util {

/*!
 * @brief Remove vectors considered null from the linear transformation matrix.
 * @param lin_trans linear transformation matrix (row-wise)
 * @param norm norms of transformed vectors
 * @param start start row to start screening
 * @param end exclusive end row to stop screening
 * @param norm_thresh threshold to decided if the vector is null
 */
template <typename value_type, typename value_type_abs>
void remove_null_vectors(subspace::Matrix<value_type>& lin_trans, std::vector<value_type_abs>& norm, const size_t start,
                         const size_t end, const value_type_abs norm_thresh) {
  for (size_t i = start, j = start; i < end; ++i) {
    if (norm[j] < norm_thresh) {
      norm.erase(begin(norm) + j);
      lin_trans.remove_row(j);
    } else
      ++j;
  }
}

/*!
 * @brief Returns true if iterator range contains integers starting with value_start and sequentially incremented by 1,
 * as in std::iota
 */
template <class ForwardIt, class EndIterator, typename Int>
bool is_iota(ForwardIt it_start, EndIterator it_end, Int value_start) {
  bool result = true;
  for (auto it = it_start; it != it_end && result; ++it, ++value_start)
    result = (*it) == value_start;
  return result;
}

} // namespace util
} // namespace itsolv
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_UTIL_H
