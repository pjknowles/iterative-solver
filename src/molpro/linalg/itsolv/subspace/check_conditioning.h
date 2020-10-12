#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_CHECK_CONDITIONING_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_CHECK_CONDITIONING_H
#include <molpro/linalg/itsolv/helper.h>
#include <molpro/linalg/itsolv/subspace/PSpace.h>
#include <molpro/linalg/itsolv/subspace/QSpace.h>
#include <molpro/linalg/itsolv/subspace/XSpaceI.h>
#include <molpro/linalg/itsolv/subspace/gram_schmidt.h>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace subspace {
namespace xspace {

/*!
 * @brief Uses Gram Schmidt procedure to transform to an orthogonal basis, removing parameters with norm less than
 * threshold from the X space.
 * @param xs subspace container
 * @param lin_trans linear transformation matrix to well conditioned basis
 * @param norm_threshold threshold for removing parameters that become too small
 * @param logger logger
 */
template <class R, class Q, class P, typename value_type, typename value_type_abs>
auto check_conditioning_gram_schmidt(XSpaceI<R, Q, P>& xspace, Matrix<value_type>& lin_trans,
                                     const value_type_abs norm_threshold, Logger& logger) {
  logger.msg("subspace::xspace::check_conditioning_gram_schmidt() ", Logger::Trace);
  auto norm = util::gram_schmidt(xspace.data[EqnData::S], lin_trans, norm_threshold);
  if (logger.data_dump)
    logger.msg("norm after Gram-Schmidt = ", begin(norm), end(norm), Logger::Info);
  auto remove_qindex = [&norm, &xspace, norm_threshold]() {
    const auto& dims = xspace.dimensions();
    auto it_begin = begin(norm) + dims.oQ;
    auto it_end = begin(norm) + dims.oQ + dims.nQ;
    auto it = std::find_if(it_begin, it_end, [&](auto el) { return el < norm_threshold; });
    return std::distance(it_begin, it);
  };
  for (auto qremove = remove_qindex(); qremove != xspace.dimensions().nQ; qremove = remove_qindex()) {
    xspace.eraseq(qremove);
    const auto ix = xspace.dimensions().oQ + qremove;
    norm.erase(begin(norm) + ix);
    lin_trans.remove_row_col(ix, ix);
    logger.msg("removed parameter index =  " + std::to_string(ix), Logger::Info);
  }
  return norm;
}

} // namespace xspace
} // namespace subspace
} // namespace itsolv
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_CHECK_CONDITIONING_H
