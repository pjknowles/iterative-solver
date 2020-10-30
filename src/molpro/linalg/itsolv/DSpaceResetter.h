#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_DSPACERESETTER_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_DSPACERESETTER_H
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/subspace/QSpace.h>
#include <molpro/linalg/itsolv/subspace/XSpaceI.h>
#include <molpro/linalg/itsolv/subspace/gram_schmidt.h>
#include <molpro/linalg/itsolv/subspace/util.h>
#include <molpro/linalg/itsolv/util.h>

namespace molpro::linalg::itsolv::detail {
//! Utility container managing whether D space resetting is in process
struct DoReset {
  DoReset() = default;
  DoReset(size_t nreset) : m_nreset(nreset) {}

  /*!
   * @brief Update status of resetting
   * @param iter current iteration counting from 0
   * @param max_size_qspace maximum size of Q space threshold stored in the solver
   * @param nD size of D space
   * @return whether resetting is on
   */
  bool update(const size_t iter, unsigned int& max_size_qspace, const size_t nD) {
    if ((iter + 1) % m_nreset == 0 && nD > 0) {
      m_value = true;
      m_save_max_size_qspace = max_size_qspace;
      max_size_qspace = std::max((unsigned int)(max_size_qspace + nD), max_size_qspace);
    }
    if (m_value && nD == 0) {
      m_value = false;
      max_size_qspace = m_save_max_size_qspace;
    }
    return m_value;
  }

  //! Whether D space is in the process of being reset
  bool value() const { return m_value; }

  //! Set the period of iterations for initiating the reset
  void set_nreset(size_t nreset) { m_nreset = nreset; }
  auto get_nreset() const { return m_nreset; }

protected:
  unsigned int m_nreset = std::numeric_limits<unsigned int>::max(); //!< reset D space every n iterations
  bool m_value = false;                                             //!< whether in the process of resetting the D space
  unsigned int m_save_max_size_qspace =
      std::numeric_limits<unsigned int>::max(); //!< stores original value of max_size_qspace before it was modified by
  //!< resetting
};

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
struct DSpaceResetter {
  DoReset do_reset;             //!< signals when resetting in in process
  std::list<Q> solution_params; //!< all current solutions that will be moved to the Q space

  template <class R, class P, typename value_type, typename value_type_abs>
  std::vector<unsigned int> run(std::vector<R>& rparams, subspace::XSpaceI<R, Q, P>& xspace,
                                const subspace::Matrix<value_type>& solutions, value_type_abs norm_thresh_solutions,
                                value_type_abs norm_thresh_null, ArrayHandlers<R, Q, P>& handlers, Logger& logger) {
    logger.msg("DSpaceResetter::run()", Logger::Trace);
    logger.msg("dimensions {nP, nQ, nD, nR} = " + std::to_string(xspace.dimensions().nP) + ", " +
                   std::to_string(xspace.dimensions().nQ) + ", " + std::to_string(xspace.dimensions().nD) + ", " +
                   std::to_string(rparams.size()),
               Logger::Trace);
    const auto nC = solutions.rows();
    auto dims = xspace.dimensions();
    if (solution_params.empty()) {
      logger.msg("constructing solutions", Logger::Debug);
      for (size_t i = 0; i < nC; ++i)
        solution_params.emplace_back(construct_zeroed_Q(rparams.first()), handlers.qr());
      auto roots = std::vector<unsigned int>(nC);
      std::iota(begin(roots), end(roots), 0);
      util::construct_solutions(wrap(begin(solution_params), end(solution_params)), roots, solutions, xspace.cparamsp(),
                                xspace.cparamsq(), xspace.cparamsd(), dims.oP, dims.oQ, dims.oD, handlers.qq(),
                                handlers.qp(), handlers.qq());
      xspace.update_dspace({}, {}, {});
    }
    const auto nR = std::min(rparams.size(), solution_params.size());
    for (size_t i = 0; i < nR; ++i) {
      handlers.rq().copy(rparams[i], solution_params.front());
      solution_params.pop_front();
    }
    const auto wparams = cwrap<R>(begin(rparams), begin(rparams) + nR);
    const auto overlap = subspace::util::overlap(wparams, xspace.cparamsq(), handlers.rq());
    for (size_t i = 0; i < nR && xspace.dimensions().nQ > 0; ++i) {
      auto* row_start = &overlap(i, xspace.dimensions().oQ);
      auto* row_end = row_start + xspace.dimensions().nQ;
      auto it_max = std::max_element(row_start, row_end);
      auto i_max = std::distance(row_start, it_max);
      xspace.eraseq(i_max);
      overlap.remove_col(xspace.dimensions().oQ + i_max);
      logger.msg("removed q index = " + std::to_string(i_max) + ", with overlap = " + std::to_string(*it_max));
    }
    // FIXME It might be prudent to double check that the subspace is stable and remove more Q vectors if it is not
    auto new_working_set = std::vector<unsigned int>(nR);
    std::iota(begin(new_working_set), end(new_working_set), 0);
    return new_working_set;
  }
};

} // namespace molpro::linalg::itsolv::detail
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_DSPACERESETTER_H
