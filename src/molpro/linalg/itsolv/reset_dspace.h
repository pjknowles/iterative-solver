#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_RESET_DSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_RESET_DSPACE_H
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/subspace/QSpace.h>
#include <molpro/linalg/itsolv/subspace/XSpaceI.h>
#include <molpro/linalg/itsolv/subspace/gram_schmidt.h>
#include <molpro/linalg/itsolv/subspace/util.h>
#include <molpro/linalg/itsolv/util.h>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace detail {

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
 * @brief Constructs solutions projected on to D space and their overlap
 * @param solutions matrix of solutions in the full subspace P+Q+D
 * @param overlap overlap matrix of the current subspace
 * @param dims dimensions of the current subspace
 * @return projection of solutions on to D space and their overlap
 */
template <typename value_type>
auto construct_projected_solutions(const subspace::Matrix<value_type>& solutions,
                                   const subspace::Matrix<value_type>& overlap,
                                   const subspace::xspace::Dimensions& dims) {
  const auto nSol = solutions.rows();
  auto overlap_DD = subspace::Matrix<value_type>({dims.nD, dims.nD});
  overlap_DD.slice() = overlap.slice({dims.oD, dims.oD}, {dims.oD + dims.nD, dims.oD + dims.nD});
  auto projected_solutions = subspace::Matrix<value_type>({nSol, dims.nD});
  projected_solutions.slice() = solutions.slice({0, dims.oD}, {nSol, dims.oD + dims.nD});
  //
  auto overlap_SS = subspace::Matrix<value_type>({nSol, nSol});
  for (size_t i = 0; i < nSol; ++i) {
    for (size_t j = 0; j <= i; ++j) {
      for (size_t k = 0; k < dims.nD; ++k) {
        for (size_t l = 0; l < dims.nD; ++l) {
          overlap_SS(i, j) += projected_solutions(i, k) * projected_solutions(j, l) * overlap_DD(k, l);
        }
      }
      overlap_SS(j, i) = overlap_SS(i, j);
    }
  }
  return std::make_tuple(projected_solutions, overlap_SS);
}

/*!
 * @brief Orthogonalise the projected solutions and transform the result to the D space
 * @param projected_solutions matrix of solutions projected on to the D space
 * @param overlap overlap matrix of projected solutions
 * @return
 */
template <typename value_type>
auto construct_orthogonalised_proj_solutions_in_D(const subspace::Matrix<value_type>& projected_solutions,
                                                  const subspace::Matrix<value_type>& overlap) {
  const auto nSol = projected_solutions.rows(), nD = projected_solutions.cols();
  auto orth_proj_sol = subspace::Matrix<value_type>{};
  auto norm = subspace::util::gram_schmidt(overlap, orth_proj_sol, 0);
  auto orth_proj_sol_D = subspace::Matrix<value_type>({nSol, nD});
  for (size_t i = 0; i < nSol; ++i)
    for (size_t j = 0; j < nD; ++j)
      for (size_t k = 0; k < nSol; ++k)
        orth_proj_sol_D(i, j) += orth_proj_sol(i, k) * projected_solutions(k, j);
  return std::make_tuple(orth_proj_sol_D, norm);
}

/*!
 * @brief Returns indices of solutions that can be added fully, and removes corresponding projected solutions
 * @param orth_proj_sol matrix of orthogonalised projected solutions
 * @param norm norms of orth_proj_sol
 * @param nR number of R parameters that can be added as full solutions
 * @param norm_thresh threshold
 * @return
 */
template <typename value_type, typename value_type_abs>
auto select_full_solutions(subspace::Matrix<value_type>& orth_proj_sol, std::vector<value_type_abs>& norm,
                           const size_t nR, const value_type_abs norm_thresh) {
  const auto nS = orth_proj_sol.rows();
  auto full_solutions = std::vector<unsigned int>{};
  for (size_t i = 0, j = 0; i < nS; ++i)
    if (norm[j] > norm_thresh && full_solutions.size() < nR) {
      full_solutions.emplace_back(i);
      norm.erase(begin(norm) + j);
      orth_proj_sol.remove_row(j);
    } else {
      ++j;
    }
  return full_solutions;
}

/*!
 * @brief Construct the rest of R space
 * @param params_remaining wrappers to the remaining R parameter containers
 * @param lin_trans remaining R space in (current R + D) subspace
 * @param rparams current R parameters
 * @param dparams D space parameters
 * @param handler_rr array handler
 * @param handler_rq array handler
 */
template <class R, class Q, typename value_type>
void construct_remaining_R(VecRef<R>& params_remaining, const subspace::Matrix<value_type>& lin_trans,
                           const CVecRef<Q>& dparams, array::ArrayHandler<R, R>& handler_rr,
                           array::ArrayHandler<R, Q>& handler_rq) {
  for (auto& param : params_remaining)
    handler_rr.fill(0, param);
  const auto nD = dparams.size();
  for (size_t i = 0; i < params_remaining.size(); ++i) {
    for (size_t j = 0; j < nD; ++j)
      handler_rq.axpy(lin_trans(i, j), dparams[j], params_remaining[i]);
  }
  for (auto& param : params_remaining) {
    auto norm = handler_rr.dot(param, param);
    norm = std::sqrt(std::abs(norm));
    handler_rr.scal(1. / norm, param);
  }
}

/*!
 * @brief Constructs the new D parameters from the remainder
 * @param lin_trans new D parameters in (R + D) subpsace
 * @param rparams R parameters
 * @param dparams D parameters
 * @param handler_rr array handler
 * @param handler_rq array handler
 * @return new D parameters and actions
 */
template <class R, class Q, typename value_type>
auto construct_new_D(const subspace::Matrix<value_type>& lin_trans, const CVecRef<Q>& dparams,
                     const CVecRef<Q>& dactions, array::ArrayHandler<Q, Q>& handler_qq,
                     array::ArrayHandler<Q, R>& handler_qr) {
  const auto nDnew = lin_trans.rows();
  const auto nD = dparams.size();
  std::vector<Q> dparams_new, dactions_new;
  {
    auto qzero = handler_qq.copy(dparams.front());
    handler_qq.fill(0, qzero);
    for (size_t i = 0; i < nDnew; ++i) {
      dparams_new.emplace_back(handler_qq.copy(qzero));
      dactions_new.emplace_back(handler_qq.copy(qzero));
    }
  }
  for (size_t i = 0; i < nDnew; ++i) {
    for (size_t j = 0; j < nD; ++j) {
      handler_qq.axpy(lin_trans(i, j), dparams.at(j), dparams_new[i]);
      handler_qq.axpy(lin_trans(i, j), dactions.at(j), dactions_new[i]);
    }
  }
  for (size_t i = 0; i < nDnew; ++i) {
    auto norm = handler_qq.dot(dparams_new[i], dparams_new[i]);
    norm = std::sqrt(std::abs(norm));
    handler_qq.scal(1. / norm, dparams_new[i]);
    handler_qq.scal(1. / norm, dactions_new[i]);
  }
  return std::make_tuple(dparams_new, dactions_new);
}

/*!
 * @brief Removes part of D space and uses it as R space so that the exact action can be calculated
 *
 * D space accumulates round-off error and has to be periodically reset. This entails using part of the subspace
 * as the new R space, and removing from D. This should be repeated until D space is empty.
 *
 * @param solver linear eigensystem solver
 * @param do_reset_dspace controls when resetting of D space is in progress
 * @param parameters parameters to propose as R space
 * @param xspace current X space container
 * @param solutions current solutions in the subspace stored as rows
 * @param norm_thresh_solutions threshold for selecting whether the full solution can be added as part of resetting
 * @param norm_thresh_null parameters with norm less than threshold are considered part of the null space
 * @param handlers array handlers
 * @param logger logger
 * @return new working set, only indicating number of parameters that need their action evaluated
 */
template <class R, class Q, class P, typename value_type, typename value_type_abs>
auto reset_dspace(LinearEigensystem<R, Q, P>& solver, std::vector<R>& parameters, subspace::XSpaceI<R, Q, P>& xspace,
                  const subspace::Matrix<value_type>& solutions, value_type_abs norm_thresh_solutions,
                  value_type_abs norm_thresh_null, ArrayHandlers<R, Q, P>& handlers, Logger& logger) {
  logger.msg("reset_dspace()", Logger::Trace);
  const auto nC = solutions.rows();
  const auto nD = xspace.dimensions().nD;
  const auto nR = std::min(parameters.size(), xspace.dimensions().nD);
  assert(nC >= nD);
  subspace::Matrix<value_type> projected_solutions, overlap;
  auto orth_proj_sol = subspace::Matrix<value_type>{};
  auto norm = std::vector<value_type_abs>{};
  std::tie(projected_solutions, overlap) =
      construct_projected_solutions(solutions, xspace.data.at(subspace::EqnData::S), xspace.dimensions());
  std::tie(orth_proj_sol, norm) = construct_orthogonalised_proj_solutions_in_D(projected_solutions, overlap);
  auto full_solutions = select_full_solutions(orth_proj_sol, norm, nR, norm_thresh_solutions);
  const auto nRfull = full_solutions.size();
  // reorder solutions so that full solutions come first.
  // This ensures that they are orthogonal to the full solutions (assuming D is orthogonal to P+Q)
  if (!util::is_iota(std::begin(full_solutions), std::end(full_solutions), 0)) {
    auto new_order = std::vector<unsigned int>(nC);
    std::iota(begin(new_order), end(new_order), 0);
    for (size_t i = 0; i < nRfull; ++i)
      std::swap(new_order[i], new_order.at(full_solutions[i]));
    std::sort(begin(new_order) + nRfull, end(new_order));
    auto solutions_full_first = subspace::Matrix<value_type>(solutions.dimensions());
    for (size_t i = 0; i < nC; ++i)
      solutions_full_first.row(i) = solutions.row(new_order[i]);
    std::tie(projected_solutions, overlap) =
        construct_projected_solutions(solutions_full_first, xspace.data.at(subspace::EqnData::S), xspace.dimensions());
    std::tie(orth_proj_sol, norm) = construct_orthogonalised_proj_solutions_in_D(projected_solutions, overlap);
    auto full_solutions_new_order = select_full_solutions(orth_proj_sol, norm, nR, norm_thresh_solutions);
    assert(util::is_iota(std::begin(full_solutions_new_order), std::end(full_solutions_new_order), 0));
    for (size_t i = 0; i < nRfull; ++i)
      full_solutions[i] = new_order.at(full_solutions_new_order[i]);
  }
  solver.solution_params(full_solutions, parameters);
  util::remove_null_vectors(orth_proj_sol, norm, 0, orth_proj_sol.rows(), norm_thresh_null);
  const auto nRremaining = std::min(nR - nRfull, orth_proj_sol.rows());
  const auto nDnew = orth_proj_sol.rows() - nRremaining;
  auto lin_trans_R = subspace::Matrix<value_type>({nRremaining, nD});
  auto lin_trans_D = subspace::Matrix<value_type>({nDnew, nD});
  lin_trans_R.slice({0, 0}, {nRremaining, nD}) = orth_proj_sol.slice({0, 0}, {nRremaining, nD});
  lin_trans_D.slice({0, 0}, {nDnew, nD}) = orth_proj_sol.slice({nRremaining, 0}, {nRremaining + nDnew, nD});
  auto params_remaining = wrap<R>(begin(parameters) + nRfull, begin(parameters) + nRfull + nRremaining);
  construct_remaining_R(params_remaining, lin_trans_R, xspace.cparamsd(), handlers.rr(), handlers.rq());
  std::vector<Q> dparams, dactions;
  std::tie(dparams, dactions) =
      construct_new_D(lin_trans_D, xspace.cparamsd(), xspace.cactionsd(), handlers.qq(), handlers.qr());
  auto lin_trans_R_component = subspace::Matrix<value_type>({nDnew, nRfull + nRremaining});
  auto wdparams = wrap(dparams);
  auto wdactions = wrap(dactions);
  xspace.update_dspace(wdparams, wdactions, lin_trans_R_component);
  const auto& working_set = solver.working_set();
  auto new_working_set =
      std::vector<unsigned int>(std::begin(working_set), std::begin(working_set) + nRfull + nRremaining);
  return new_working_set;
}

} // namespace detail
} // namespace itsolv
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_RESET_DSPACE_H
