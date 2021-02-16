#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_UTIL_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_UTIL_H
#include <molpro/linalg/itsolv/ArrayHandlers.h>
#include <molpro/linalg/itsolv/subspace/Matrix.h>
#include <molpro/linalg/itsolv/wrap.h>
#include <molpro/linalg/itsolv/Statistics.h>

#include <locale>
#include <string>
#include <map>

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
  std::sort(std::begin(indices), std::end(indices), std::greater<int>());
  for (auto i : indices)
    params.erase(begin(params) + i);
}

//! Wraps some useful string manipulation functions
class StringFacet {
public:
  std::string toupper(std::string in) const;
  std::string tolower(std::string in) const;
  bool tobool(const std::string& in) const;
  static void crop_space(std::string& path);
  static std::map<std::string,std::string> parse_keyval_string(std::string s);

private:
  const std::ctype<char>& facet = std::use_facet<std::ctype<char>>(std::locale());
};

template <class R, class Q, class P>
void accumulate_handler_calls(std::shared_ptr<Statistics> stats, std::shared_ptr<ArrayHandlers<R,Q,P>> handlers) {
  stats->r_scal += handlers->rr().counter().scal;
  stats->r_scal += handlers->rq().counter().scal;
  stats->r_scal += handlers->rp().counter().scal;
  stats->rr_dot += handlers->rr().counter().dot;
  stats->rr_axpy += handlers->rr().counter().axpy;
  stats->rr_copy += handlers->rr().counter().copy;
  stats->rr_gemm_inner += handlers->rr().counter().gemm_inner;
  stats->rr_gemm_outer += handlers->rr().counter().gemm_outer;
  stats->q_scal += handlers->qq().counter().scal;
  stats->q_scal += handlers->qr().counter().scal;
  stats->q_scal += handlers->qp().counter().scal;
  stats->qq_dot += handlers->qq().counter().dot;
  stats->qq_axpy += handlers->qq().counter().axpy;
  stats->qq_copy += handlers->qq().counter().copy;
  stats->qq_gemm_inner += handlers->qq().counter().gemm_inner;
  stats->qq_gemm_outer += handlers->qq().counter().gemm_outer;
  stats->qr_dot += handlers->qr().counter().dot;
  stats->qr_axpy += handlers->qr().counter().axpy;
  stats->qr_copy += handlers->qr().counter().copy;
  stats->qr_gemm_inner += handlers->qr().counter().gemm_inner;
  stats->qr_gemm_outer += handlers->qr().counter().gemm_outer;
  stats->rq_dot += handlers->rq().counter().dot;
  stats->rq_axpy += handlers->rq().counter().axpy;
  stats->rq_copy += handlers->rq().counter().copy;
  stats->rq_gemm_inner += handlers->rq().counter().gemm_inner;
  stats->rq_gemm_outer += handlers->rq().counter().gemm_outer;
  stats->rp_dot += handlers->rp().counter().dot;
  stats->rp_axpy += handlers->rp().counter().axpy;
  stats->rp_gemm_inner += handlers->rp().counter().gemm_inner;
  stats->rp_gemm_outer += handlers->rp().counter().gemm_outer;
  stats->qp_dot += handlers->qp().counter().dot;
  stats->qp_axpy += handlers->qp().counter().axpy;
  stats->qp_gemm_inner += handlers->qp().counter().gemm_inner;
  stats->qp_gemm_outer += handlers->qp().counter().gemm_outer;
};

} // namespace molpro::linalg::itsolv::util
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_UTIL_H
