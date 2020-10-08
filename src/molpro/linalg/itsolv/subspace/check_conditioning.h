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
 * @brief Generates linear transformation of X space to an orthogonal subspace and removes any redundant Q vectors.
 * @param xs X space container
 * @param rs R space container
 * @param qs Q space container
 * @param ps P space container
 * @param cs C space container
 * @param lin_trans  Linear transformation matrix to an orthogonal subspace
 * @param norm_threshold
 * @param logger
 */
template <class R, class P, class Q, class ST>
void check_conditioning_gram_schmidt(XSpaceI<R, Q, P>& xs, QSpace<R, Q, P>& qs, PSpace<R, P>& ps,
                                     CSpace<R, Q, P, ST> cs, Matrix<ST>& lin_trans, double norm_threshold,
                                     Logger& logger) {
  logger.msg("xspace::check_conditioning_gram_schmidt", Logger::Trace);
  bool stable = false;
  auto candidates = std::vector<size_t>{qs.size()};
  std::iota(begin(candidates), end(candidates), size_t{xs.dimensions().oQ});
  while (!stable && !candidates.empty()) {
    const auto& dim = xs.dimensions();
    const auto& s = xs.data[EqnData::S];
    auto norm = util::gram_schmidt(s, lin_trans);
    if (logger.data_dump)
      logger.msg("norm after Gram-Schmidt = ", begin(norm), end(norm), Logger::Info);
    auto imin = std::find_if(begin(candidates), end(candidates),
                             [&norm, norm_threshold](auto c) { return norm[c] < norm_threshold; });
    stable = (imin == candidates.end());
    if (!stable) {
      logger.msg("removing candidate from q space i =" + std::to_string(*imin) +
                     ", norm = " + Logger::scientific(norm[*imin]),
                 Logger::Debug);
      qs.erase(*imin);
      xs.build_subspace(qs, ps, cs);
      candidates.resize(qs.size());
    }
  }
}

template <class R, class Q, class P, typename value_type, typename value_type_abs>
void check_conditioning_gram_schmidt(XSpaceI<R, Q, P>& xs, Matrix<value_type>& lin_trans,
                                     const value_type_abs norm_threshold, Logger& logger) {}

} // namespace xspace
} // namespace subspace
} // namespace itsolv
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_CHECK_CONDITIONING_H
