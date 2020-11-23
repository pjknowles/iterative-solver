#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_DSPACERESETTER_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_DSPACERESETTER_H
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/propose_rspace.h>
#include <molpro/linalg/itsolv/subspace/Dimensions.h>
#include <molpro/linalg/itsolv/subspace/QSpace.h>
#include <molpro/linalg/itsolv/subspace/XSpaceI.h>
#include <molpro/linalg/itsolv/subspace/gram_schmidt.h>
#include <molpro/linalg/itsolv/subspace/util.h>
#include <molpro/linalg/itsolv/util.h>

namespace molpro::linalg::itsolv::detail {
//! Removes Q parameters that have smallest contribution to any solution until Q space size is within limit
template <class R, class Q, class P, typename value_type>
void resize_qspace(subspace::XSpaceI<R, Q, P>& xspace, const subspace::Matrix<value_type>& solutions,
                   int m_max_Qsize_after_reset, Logger& logger) {
  logger.msg("resize_qspace()", Logger::Trace);
  auto q_delete = limit_qspace_size(xspace.dimensions(), m_max_Qsize_after_reset, solutions, logger);
  logger.msg("delete Q parameter indices = ", q_delete.begin(), q_delete.end(), Logger::Debug);
  std::sort(begin(q_delete), end(q_delete), std::greater());
  for (auto iq : q_delete)
    xspace.eraseq(iq);
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
  int m_nreset = std::numeric_limits<int>::max();                //!< reset D space every n iterations
  int m_max_Qsize_after_reset = std::numeric_limits<int>::max(); //!< maximum size of Q space after reset
  std::list<Q> solution_params; //!< all current solutions that will be moved to the Q space

public:
  DSpaceResetter() = default;
  DSpaceResetter(int nreset, int max_Qsize) : m_nreset{nreset}, m_max_Qsize_after_reset{max_Qsize} {}

  //! Whether reset operation should be run
  bool do_reset(size_t iter, const subspace::Dimensions& dims) {
    return ((iter + 1) % m_nreset == 0 && dims.nD > 0) || !solution_params.empty();
  }

  void set_nreset(size_t i) {
    assert(i >= 0);
    m_nreset = i;
  }
  auto get_nreset() const { return m_nreset; }
  void set_max_Qsize(size_t i) {
    assert(i >= 0);
    m_max_Qsize_after_reset = i;
  }
  auto get_max_Qsize() const { return m_max_Qsize_after_reset; }

  //! Run the reset operation
  template <class R, class P, typename value_type>
  std::vector<int> run(const VecRef<R>& rparams, subspace::XSpaceI<R, Q, P>& xspace,
                       const subspace::Matrix<value_type>& solutions, ArrayHandlers<R, Q, P>& handlers,
                       Logger& logger) {
    logger.msg("DSpaceResetter::run()", Logger::Trace);
    logger.msg("dimensions {nP, nQ, nD, nR} = " + std::to_string(xspace.dimensions().nP) + ", " +
                   std::to_string(xspace.dimensions().nQ) + ", " + std::to_string(xspace.dimensions().nD) + ", " +
                   std::to_string(rparams.size()),
               Logger::Trace);
    const auto nC = solutions.rows();
    auto dims = xspace.dimensions();
    if (solution_params.empty() && !rparams.empty()) {
      logger.msg("constructing solutions", Logger::Debug);
      for (size_t i = 0; i < nC; ++i)
        solution_params.emplace_back(util::construct_zeroed_copy(rparams.front().get(), handlers.qr()));
      auto roots = std::vector<int>(nC);
      std::iota(begin(roots), end(roots), 0);
      util::construct_solutions(wrap<Q>(begin(solution_params), end(solution_params)), roots, solutions,
                                xspace.cparamsp(), xspace.cparamsq(), xspace.cparamsd(), dims.oP, dims.oQ, dims.oD,
                                handlers.qq(), handlers.qp(), handlers.qq());
      VecRef<Q> null_params, null_actions;
      xspace.update_dspace(null_params, null_actions);
    }
    const auto nR = std::min(rparams.size(), solution_params.size());
    for (size_t i = 0; i < nR; ++i) {
      handlers.rq().copy(rparams[i], solution_params.front());
      solution_params.pop_front();
    }
    const auto wparams = cwrap<R>(begin(rparams), begin(rparams) + nR);
    auto overlap = subspace::util::overlap(wparams, xspace.cparamsq(), handlers.rq());
    for (size_t i = 0; i < nR && xspace.dimensions().nQ > 0; ++i) {
      auto* row_start = &overlap(i, xspace.dimensions().oQ);
      auto* row_end = row_start + xspace.dimensions().nQ;
      auto it_max = std::max_element(row_start, row_end);
      auto i_max = std::distance(row_start, it_max);
      xspace.eraseq(i_max);
      overlap.remove_col(xspace.dimensions().oQ + i_max);
      logger.msg("removed q index = " + std::to_string(i_max) + ", with overlap = " + std::to_string(*it_max),
                 Logger::Debug);
    }
    resize_qspace(xspace, solutions, m_max_Qsize_after_reset, logger);
    auto new_working_set = std::vector<int>(nR);
    std::iota(begin(new_working_set), end(new_working_set), 0);
    return new_working_set;
  }
};

} // namespace molpro::linalg::itsolv::detail
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_DSPACERESETTER_H
