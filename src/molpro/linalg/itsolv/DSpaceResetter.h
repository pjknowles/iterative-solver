#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_DSPACERESETTER_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_DSPACERESETTER_H
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/subspace/QSpace.h>
#include <molpro/linalg/itsolv/subspace/XSpaceI.h>
#include <molpro/linalg/itsolv/subspace/gram_schmidt.h>
#include <molpro/linalg/itsolv/subspace/util.h>
#include <molpro/linalg/itsolv/util.h>

namespace molpro::linalg::itsolv::detail {
/*!
 * @brief Resets D space constructing full solutions as the new working set, removing instabilities from Q space, and
 * clearing D space
 *
 * Resetting can happen over multiple iterations if size of the working set is less than the number of solutions.
 * In such cases, D space is still cleared at the start, so only the first batch of solutions will be reproduced at
 * that iteration with the other possibly showing larger errors than before. By the end of resetting, all previous
 * solutions will be in the Q space and previous energies and residuals will be reproduced.
 *
 * @tparam Q
 */
template <class Q>
class DSpaceResetter {
protected:
  unsigned int m_nreset = std::numeric_limits<unsigned int>::max(); //!< reset D space every n iterations
  std::list<Q> solution_params; //!< all current solutions that will be moved to the Q space

public:
  //! Whether reset operation should be run
  bool do_reset(size_t iter, const subspace::xspace::Dimensions& dims) {
    return ((iter + 1) % m_nreset == 0 && dims.nD > 0) || !solution_params.empty();
  }

  void set_nreset(size_t i) { m_nreset = i; }
  auto get_nreset() const { return m_nreset; }

  //! Run the reset operation
  template <class R, class P, typename value_type>
  std::vector<unsigned int> run(std::vector<R>& rparams, subspace::XSpaceI<R, Q, P>& xspace,
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
        solution_params.emplace_back(util::construct_zeroed_copy(rparams.front(), handlers.qr()));
      auto roots = std::vector<unsigned int>(nC);
      std::iota(begin(roots), end(roots), 0);
      util::construct_solutions(wrap<Q>(begin(solution_params), end(solution_params)), roots, solutions,
                                xspace.cparamsp(), xspace.cparamsq(), xspace.cparamsd(), dims.oP, dims.oQ, dims.oD,
                                handlers.qq(), handlers.qp(), handlers.qq());
      VecRef<Q> null_params, null_actions;
      xspace.update_dspace(null_params, null_actions, {});
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
    // FIXME It might be prudent to double check that the subspace is stable and remove more Q vectors if it is not
    auto new_working_set = std::vector<unsigned int>(nR);
    std::iota(begin(new_working_set), end(new_working_set), 0);
    return new_working_set;
  }
};

} // namespace molpro::linalg::itsolv::detail
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_DSPACERESETTER_H
