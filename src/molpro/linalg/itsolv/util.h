#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_UTIL_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_UTIL_H
#include <molpro/linalg/itsolv/ArrayHandlers.h>
#include <molpro/linalg/itsolv/subspace/Matrix.h>
#include <molpro/linalg/itsolv/wrap.h>

#include <locale>
#include <string>

namespace molpro::linalg::itsolv::util {

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

//! Copy R parameter to Q  and zero it. This ensures same size and distribution of the copy.
template <class R, class Q>
Q construct_zeroed_copy(const R& param, array::ArrayHandler<Q, R>& handler) {
  auto q = handler.copy(param);
  handler.fill(0, q);
  return q;
}

/*!
 * @brief Construct solutions
 * @param params parameters for storing solutions
 * @param roots roots to construct
 * @param solutions solution matrix in the subspace with solutions in rows
 * @param pparams P space parameters
 * @param qparams Q space parameters
 * @param dparams D space parameters
 * @param oP offset to P space in subspace
 * @param oQ offset to Q space in subspace
 * @param oD offset to D space in subspace
 * @param handler_rr array handler
 * @param handler_rp array handler
 * @param handler_rq array handler
 */
template <class R, class Q, class P>
void construct_solutions(const VecRef<R>& params, const std::vector<int>& roots,
                         const subspace::Matrix<double>& solutions, const CVecRef<P>& pparams,
                         const CVecRef<Q>& qparams, const CVecRef<Q>& dparams, size_t oP, size_t oQ, size_t oD,
                         array::ArrayHandler<R, R>& handler_rr, array::ArrayHandler<R, P>& handler_rp,
                         array::ArrayHandler<R, Q>& handler_rq) {
  assert(params.size() >= roots.size());
  for (size_t i = 0; i < roots.size(); ++i)
    handler_rr.fill(0, params.at(i));
  for (size_t i = 0; i < roots.size(); ++i) {
    auto root = roots[i];
    for (size_t j = 0; j < pparams.size(); ++j)
      handler_rp.axpy(solutions(root, oP + j), pparams.at(j), params.at(i));
    for (size_t j = 0; j < qparams.size(); ++j)
      handler_rq.axpy(solutions(root, oQ + j), qparams.at(j), params.at(i));
    for (size_t j = 0; j < dparams.size(); ++j)
      handler_rq.axpy(solutions(root, oD + j), dparams.at(j), params.at(i));
  }
}

/*!
 * @brief Removes parameters
 * @param indices indices of parameters to remove
 * @param params list of parameters
 */
template <class Container>
void delete_parameters(std::vector<int> indices, Container& params) {
  std::sort(std::begin(indices), std::end(indices), std::greater());
  for (auto i : indices)
    params.erase(begin(params) + i);
}

//! Wraps some useful string manipulation functions
class StringFacet {
public:
  std::string toupper(std::string in);
  std::string tolower(std::string in);
  bool tobool(const std::string& in);
  static void crop_space(std::string& path);

private:
  const std::ctype<char>& facet = std::use_facet<std::ctype<char>>(std::locale());
};

} // namespace molpro::linalg::itsolv::util
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_UTIL_H
