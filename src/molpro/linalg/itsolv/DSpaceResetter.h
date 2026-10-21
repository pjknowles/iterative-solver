#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_DSPACERESETTER_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_DSPACERESETTER_H
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/propose_rspace.h>
#include <molpro/linalg/itsolv/qspace_options.h>
#include <molpro/linalg/itsolv/subspace/Dimensions.h>
#include <molpro/linalg/itsolv/subspace/IXSpace.h>
#include <molpro/linalg/itsolv/subspace/QSpace.h>
#include <molpro/linalg/itsolv/subspace/gram_schmidt.h>
#include <molpro/linalg/itsolv/subspace/util.h>
#include <molpro/linalg/itsolv/util.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <limits>

namespace molpro::linalg::itsolv::detail {
//! Removes Q parameters that have smallest contribution to any solution until Q space size is within limit
template <class R, class Q, class P, typename value_type>
void resize_qspace(subspace::IXSpace<R, Q, P>& xspace, const subspace::Matrix<value_type>& solutions,
                   std::size_t max_size, Logger& logger) {
  logger.trace("resize_qspace()");
  auto q_delete = limit_qspace_size(xspace.dimensions(), QSpaceOptions{.max_size = max_size}, solutions, logger);
  logger.debug("delete Q parameter indices = ", q_delete);
  erase_qspace(xspace, q_delete);
}

//! Erases Q parameters with the given indices together with the corresponding columns of a solutions matrix
template <class R, class Q, class P, typename value_type, typename Index>
void erase_qspace(subspace::IXSpace<R, Q, P>& xspace, subspace::Matrix<value_type>& solutions,
                  std::vector<Index> indices) {
  assert(solutions.cols() == xspace.dimensions().nX());
  const auto oQ = xspace.dimensions().oQ();
  std::sort(begin(indices), end(indices), std::greater<Index>());
  for (auto i : indices) {
    xspace.eraseq(i);
    solutions.remove_col(oQ + i);
  }
}

/*!
 * @brief Returns list of Q indices with maximum overlap to R parameters, sorted in descending order
 * @param rparams R parameters
 * @param qparams Q parameters to calculate maximum overlap of with R
 * @param handler array handler
 * @param logger logger
 */
template <class R, class Q>
auto max_overlap_with_R(const CVecRef<R>& rparams, const CVecRef<Q>& qparams, array::ArrayHandler<R, Q>& handler,
                        Logger& logger) {
  logger.trace("max_overlap_with_R()");
  const auto nR = rparams.size();
  auto overlap = subspace::util::overlap(rparams, qparams, handler);
  auto q_indices = std::vector<int>(qparams.size());
  std::iota(std::begin(q_indices), std::end(q_indices), 0);
  auto q_max_overlap = std::vector<int>{};
  for (size_t i = 0; i < nR && !q_indices.empty(); ++i) {
    auto ov = std::vector<typename array::ArrayHandler<R, Q>::value_type_abs>{};
    for (auto j : q_indices)
      ov.push_back(std::abs(overlap(i, j)));
    auto it_max = std::max_element(ov.begin(), ov.end());
    auto i_max = std::distance(ov.begin(), it_max);
    logger.debug("removed q index = " + std::to_string(q_indices[i_max]) +
                 ", with overlap = " + std::to_string(*it_max));
    q_max_overlap.push_back(q_indices[i_max]);
    q_indices.erase(q_indices.begin() + i_max);
  }
  std::sort(std::begin(q_max_overlap), std::end(q_max_overlap), std::greater<int>());
  return q_max_overlap;
}

/*!
 * @brief Resets D space constructing full solutions as the new working set, removing instabilities from Q space, and
 * clearing D space
 *
 * Resetting can happen over multiple iterations if size of the working set is less than the number of solutions.
 * In such cases, D space is still cleared at the start, so only the first batch of solutions will be reproduced at
 * that iteration with the other possibly showing larger errors than before. By the end of resetting, all previous
 * solutions will be in the Q space and previous energies and residuals will be reproduced.
 *
 * After resetting, the size of Q space can be reduced. The maximum size can be controlled via set_max_Qsize().
 *
 * @tparam Q
 */
template <class Q>
class DSpaceResetter {
protected:
  std::size_t m_nreset = 0; //!< reset D space every n iterations, 0 for never
  std::size_t m_max_Qsize_after_reset =
      std::numeric_limits<std::size_t>::max(); //!< maximum size of Q space after reset
  std::list<Q> solution_params;                //!< all current solutions that will be moved to the Q space

public:
  DSpaceResetter() = default;
  DSpaceResetter(std::size_t nreset, std::size_t max_Qsize) : m_nreset{nreset}, m_max_Qsize_after_reset{max_Qsize} {}

  //! Whether reset operation should be run
  bool do_reset(size_t iter, const subspace::Dimensions& dims) {
    return (m_nreset != 0 && (iter + 1) % m_nreset == 0 && dims.nD > 0) || !solution_params.empty();
  }

  void set_nreset(size_t i) { m_nreset = i; }
  auto get_nreset() const { return m_nreset; }
  void set_max_Qsize(size_t i) { m_max_Qsize_after_reset = i; }
  auto get_max_Qsize() const { return m_max_Qsize_after_reset; }

  //! Run the reset operation
  template <class R, class P, typename value_type, typename value_type_abs>
  std::vector<int> run(const VecRef<R>& rparams, subspace::IXSpace<R, Q, P>& xspace,
                       const subspace::Matrix<value_type>& solutions, const value_type_abs norm_thresh,
                       const value_type_abs svd_thresh, ArrayHandlers<R, Q, P>& handlers, Logger& logger) {
    logger.trace("DSpaceResetter::run()");
    logger.trace("dimensions {nP, nQ, nD, nR} = ", xspace.dimensions().nP, xspace.dimensions().nQ,
                 xspace.dimensions().nD, rparams.size());
    // solutions with columns kept in step with the changes to xspace below
    auto solutions_xspace = solutions;
    if (solution_params.empty() && !rparams.empty()) {
      logger.debug("constructing solutions");
      const auto dims = xspace.dimensions();
      const auto overlap = xspace.data.at(subspace::EqnData::S);
      auto q_indices = std::vector<std::size_t>(xspace.dimensions().nQ);
      std::iota(q_indices.begin(), q_indices.end(), 0);
      const auto solutions_proj =
          dspace::stable_projected_solutions(solutions, overlap, dims, q_indices, norm_thresh, svd_thresh, logger);
      const auto nC = solutions_proj.rows();
      for (size_t i = 0; i < nC; ++i)
        solution_params.emplace_back(util::construct_zeroed_copy(rparams.front().get(), handlers.qr()));
      auto roots = std::vector<int>(nC);
      std::iota(begin(roots), end(roots), 0);
      util::construct_solutions(wrap(begin(solution_params), end(solution_params)), roots, solutions_proj, {},
                                xspace.cparamsq(), xspace.cparamsd(), 0, 0, dims.nQ, handlers.qq(), handlers.qp(),
                                handlers.qq());
      assert(solutions_xspace.cols() == dims.nX());
      for (size_t i = 0; i < dims.nD; ++i)
        solutions_xspace.remove_col(dims.oD());
      VecRef<Q> null_params, null_actions;
      xspace.update_dspace(null_params, null_actions);
    }
    const auto nR = std::min(rparams.size(), solution_params.size());
    for (size_t i = 0; i < nR; ++i) {
      handlers.rq().copy(rparams[i], solution_params.front());
      solution_params.pop_front();
    }
    const auto wparams = cwrap(begin(rparams), begin(rparams) + nR);
    auto q_delete = max_overlap_with_R(wparams, xspace.cparamsq(), handlers.rq(), logger);
    erase_qspace(xspace, solutions_xspace, q_delete);
    assert(solutions_xspace.cols() == xspace.dimensions().nX());
    if (xspace.dimensions().nQ + nR > m_max_Qsize_after_reset)
      resize_qspace(xspace, solutions_xspace, m_max_Qsize_after_reset > nR ? m_max_Qsize_after_reset - nR : 0, logger);
    auto new_working_set = std::vector<int>(nR);
    std::iota(begin(new_working_set), end(new_working_set), 0);
    return new_working_set;
  }
};

} // namespace molpro::linalg::itsolv::detail
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_DSPACERESETTER_H
