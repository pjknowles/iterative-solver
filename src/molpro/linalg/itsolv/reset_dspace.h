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
  /*!
   * @brief Update status of resetting
   * @param iter current iteration
   * @param max_size_qspace maximum size of Q space threshold stored in the solver
   * @param nD size of D space
   * @return whether resetting is on
   */
  bool update(const size_t iter, unsigned int& max_size_qspace, const size_t nD) {
    if ((iter + 1) % n_reset_D == 0) {
      value = true;
      m_save_max_size_qspace = max_size_qspace;
      max_size_qspace = std::max((unsigned int)(max_size_qspace + nD), max_size_qspace);
    }
    if (value && nD == 0) {
      value = false;
      max_size_qspace = m_save_max_size_qspace;
    }
    return value;
  }
  unsigned int n_reset_D = std::numeric_limits<unsigned int>::max(); //!< reset D space every n iterations
  bool value = false; //!< whether in the process of resetting the D space
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
 * @brief Constructs overlap matrix for R+D subspace
 * @param rparams R space parameters
 * @param dparams D space parameters
 * @param overlap_DD overlap matrix for D subspace
 * @param handler_rr array handler
 * @param handler_rq array handler
 */
template <class R, class Q, typename value_type>
auto construct_overlap_RD(const CVecRef<R>& rparams, const CVecRef<Q>& dparams,
                          subspace::Matrix<value_type>& overlap_DD, array::ArrayHandler<R, R>& handler_rr,
                          array::ArrayHandler<R, Q>& handler_rq) {
  const auto nR = rparams.size();
  const auto nD = dparams.size();
  auto overlap = subspace::Matrix<value_type>({nR + nD, nR + nD});
  overlap.slice({0, 0}, {nR, nR}) = subspace::util::overlap(rparams, handler_rr);
  overlap.slice({0, nR}, {nR, nR + nD}) = subspace::util::overlap(rparams, dparams, handler_rq);
  overlap.slice({nR, nR}, {nR + nD, nR + nD}) = overlap_DD;
  subspace::transpose_copy(overlap.slice({nR, 0}, {nR + nD, nR}), overlap.slice({0, nR}, {nR, nR + nD}));
  return overlap;
}

/*!
 * @brief Propose transformation to the subspace of D that is orthogonal to current R
 * @param overlap_RD overlap matrix of R+D subspace
 * @param nR size of R space
 * @param norm_thresh_null threshold indicating null space
 */
template <typename value_type, typename value_type_abs>
auto propose_remaining_D_space(const subspace::Matrix<value_type>& overlap_RD, const size_t nR,
                               const value_type_abs norm_thresh_null) {
  assert(overlap_RD.rows() > nR);
  const auto nD = overlap_RD.rows() - nR;
  auto lin_trans = subspace::Matrix<value_type>(overlap_RD.dimensions());
  auto norm = subspace::util::gram_schmidt(overlap_RD, lin_trans, norm_thresh_null);
  auto wnorm = wrap(norm);
  std::sort(std::begin(wnorm) + nR, std::end(wnorm), std::greater<value_type_abs>{});
  auto order = find_ref(wnorm, begin(norm), end(norm));
  auto norm_ordered = std::vector<value_type_abs>();
  std::copy(std::begin(wnorm) + nR, std::end(wnorm), std::back_inserter(norm_ordered));
  auto lin_trans_ordered = subspace::Matrix<value_type>({nD, nR + nD});
  for (size_t i = 0; i < nD; ++i) {
    const auto ii = order.at(nR + i);
    lin_trans_ordered.row(i) = lin_trans.row(ii);
  }
  util::remove_null_vectors(lin_trans_ordered, norm_ordered, 0, lin_trans_ordered.rows(), norm_thresh_null);
  for (size_t i = 0; i < nD; ++i) {
    lin_trans_ordered.row(i).scal(1. / norm_ordered[i]);
  }
  return lin_trans_ordered;
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
                           const CVecRef<R>& rparams, const CVecRef<Q>& dparams, array::ArrayHandler<R, R>& handler_rr,
                           array::ArrayHandler<R, Q>& handler_rq) {
  for (auto& param : params_remaining)
    handler_rr.fill(0, param);
  const auto nR = rparams.size();
  const auto nD = dparams.size();
  for (size_t i = 0; i < params_remaining.size(); ++i) {
    for (size_t j = 0; j < nR; ++j)
      handler_rr.axpy(lin_trans(i, j), rparams[j], params_remaining[i]);
    for (size_t j = 0; j < nD; ++j)
      handler_rq.axpy(lin_trans(i, nR + j), dparams[j], params_remaining[i]);
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
auto construct_new_D(const subspace::Matrix<value_type>& lin_trans, const CVecRef<R>& rparams,
                     const CVecRef<Q>& dparams, const CVecRef<Q>& dactions, array::ArrayHandler<Q, Q>& handler_qq,
                     array::ArrayHandler<Q, R>& handler_qr) {
  const auto nDnew = lin_trans.rows();
  const auto nR = rparams.size(), nD = dparams.size();
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
    for (size_t j = 0; j < nR; ++j) {
      handler_qr.axpy(lin_trans(i, j), rparams.at(j), dparams_new[i]);
    }
    for (size_t j = 0; j < nD; ++j) {
      handler_qq.axpy(lin_trans(i, nR + j), dparams.at(j), dparams_new[i]);
      handler_qq.axpy(lin_trans(i, nR + j), dactions.at(j), dactions_new[i]);
    }
  }
  return std::make_tuple(dparams_new, dactions_new);
}

/*!
 * @brief Construct overlap matrix of P+Q+Solutions subspace
 * @param solutions matrix in the current subspace with solutions as rows
 * @param overlap_xx overlap matrix of current subspace (P+Q+D), which spans the solutions
 * @param dims dimensions of current subspace
 * @return
 */
template <typename value_type>
auto construct_overlap_with_solutions(const subspace::Matrix<value_type>& solutions,
                                      const subspace::Matrix<value_type>& overlap_xx,
                                      const subspace::xspace::Dimensions& dims) {
  const auto nSol = solutions.rows();
  const auto nPQ = dims.nP + dims.nQ;
  const auto oSol = nPQ;
  auto overlap = overlap_xx;
  overlap.resize({nPQ + nSol, nPQ + nSol});
  for (size_t i = 0; i < nSol; ++i) {
    for (size_t j = 0; j < nPQ; ++j) {
      value_type v = 0;
      for (size_t k = 0; k < dims.nX; ++k) {
        v += solutions(i, k) * overlap_xx(k, j);
      }
      overlap(oSol + i, j) = overlap(j, oSol + i) = v;
    }
  }
  for (size_t i = 0; i < nSol; ++i) {
    for (size_t j = 0; j <= i; ++j) {
      value_type v = 0;
      for (size_t k = 0; k < dims.nX; ++k) {
        for (size_t l = 0; l < dims.nX; ++l) {
          v += solutions(i, k) * solutions(j, l) * overlap_xx(k, l);
        }
      }
      overlap(oSol + i, oSol + j) = overlap(oSol + j, oSol + i) = v;
    }
  }
  return overlap;
}

/*!
 * @brief Given linear transformation to R parameters in subspace P+Q+Solutions, transform it to P+Q+D subspace
 *
 * x_i &= \sum_{j \in P+Q} L_ij v_j + \sum_{j \in C} T_ij z_j \\
 *     &= \sum_{j \in P+Q} L_ij v_j + \sum_{j \in C} \sum_{k \in P+Q+D} T_ij S_jk v_k \\
 *     &= \sum_{j \in P+Q} L_ij v_j + \sum_{k \in P+Q+D} D_ik v_k \\
 * D_ik &= \sum_{j \in C} T_ij S_jk
 *
 * @param lin_trans transformation to R parameters in subspace P+Q+Solutions
 * @param solutions solution matrix in current subspace P+Q+D
 * @param dims dimensions of current subspace
 * @return linear transformation to R parameters in P+Q+D subspace
 */
template <typename value_type>
auto transform_PQSol_to_PQD_subspace(const subspace::Matrix<value_type>& lin_trans,
                                     const subspace::Matrix<value_type>& solutions,
                                     const subspace::xspace::Dimensions& dims) {
  const auto nR = lin_trans.rows();
  const auto nP = dims.nP, nQ = dims.nQ, nD = dims.nD, nX = nP + nQ + nD;
  auto trans_pqd = subspace::Matrix<value_type>({nR, nP + nQ + nD});
  trans_pqd.slice({0, 0}, {nR, nP + nQ}) = lin_trans.slice({0, 0}, {nR, nP + nQ});
  for (size_t i = 0; i < nR; ++i)
    for (size_t j = 0; j < nX; ++j)
      for (size_t k = 0; k < nR; ++k)
        trans_pqd(i, j) += lin_trans(i, k) * solutions(k, j);
  return trans_pqd;
}

/*!
 * @brief Given solutions orthogonalised against P+Q subspace proposes R parameters  and new D parameters.
 *
 * The number of solutions could be greater than size of D space (nSol >= nD), and nR <= nD.
 * Thus, orthogonalised solutions are considered in order of descending norm. By construction, there is at least
 * nR + nDnew parameters that are not null.
 *
 * Parameters for R space are full solutions if the orthogonalised solution has norm greater than norm_thresh.
 *
 * @param lin_trans linear transformation to solutions orthogonalised against P+Q subspace
 * @param norm norm of proposed transformed parameters
 * @param solutions matrix of solutions in the current P+Q+D subspace
 * @param nR number of R parameters that will be proposed
 * @param nDnew number of new D parameters
 * @param norm_thresh threshold for selecting which solutions can be added fully
 */
template <typename value_type, typename value_type_abs>
auto propose_R_and_D_params(subspace::Matrix<value_type>& lin_trans, std::vector<value_type_abs>& norm,
                            const subspace::Matrix<value_type>& solutions, const size_t nR, const size_t nDnew,
                            value_type_abs norm_thresh) {
  assert(lin_trans.rows() == norm.size());
  const auto nX = lin_trans.cols();
  auto wnorm = wrap(norm);
  auto descending_norm = wnorm;
  std::sort(std::begin(descending_norm), std::end(descending_norm), std::greater<value_type_abs>{});
  auto order = find_ref(descending_norm, std::begin(norm), std::end(norm));
  auto lin_trans_R = subspace::Matrix<value_type>({nR, nX});
  for (size_t i = 0; i < nR; ++i) {
    const auto ii = order[i];
    if (norm[ii] > norm_thresh) {
      lin_trans_R.row(i) = solutions.row(ii);
    } else {
      lin_trans_R.row(i) = lin_trans.row(ii);
      lin_trans_R.row(i).scal(1. / norm[ii]);
    }
  }
  auto lin_trans_D = subspace::Matrix<value_type>({nDnew, nX});
  for (size_t i = 0; i < nDnew; ++i) {
    const auto ii = order[nR + i];
    lin_trans_D.row(i) = lin_trans.row(ii);
    lin_trans_D.row(i).scal(1. / norm[ii]);
  }
  return std::make_tuple(lin_trans_R, lin_trans_D);
}

/*!
 * @brief Constructs the proposed R and D space parameters from the current subspace (P+Q+D)
 * @param parameters containers for R parameters
 * @param lin_trans_R proposed R parameters in the current subspace
 * @param lin_trans_D proposed D parameters in the current subspace
 * @param pparams P space parameters
 * @param pactions P space actions
 * @param qparams Q space parameters
 * @param qactions Q space actions
 * @param dparams D space parameters
 * @param dactions D space actions
 * @param handlers array handlers
 * @return new D space parameters and actions
 */
template <class R, class Q, class P, typename value_type>
auto construct_R_and_D_params(std::vector<R>& parameters, const subspace::Matrix<value_type>& lin_trans_R,
                              const subspace::Matrix<value_type>& lin_trans_D, const CVecRef<P>& pparams,
                              const CVecRef<P>& pactions, const CVecRef<Q>& qparams, const CVecRef<Q>& qactions,
                              const CVecRef<Q>& dparams, const CVecRef<Q>& dactions, ArrayHandlers<R, Q, P>& handlers) {
  const auto nR = lin_trans_R.rows();
  const auto nDnew = lin_trans_D.rows();
  const auto nP = pparams.size(), nQ = qparams.size(), nD = dparams.size();
  std::vector<Q> dparams_new, dactions_new;
  {
    auto qzero = handlers.qq().copy(qparams.front());
    handlers.qq().fill(0, qzero);
    for (size_t i = 0; i < nDnew; ++i) {
      dparams_new.emplace_back(handlers.qq().copy(qzero));
      dactions_new.emplace_back(handlers.qq().copy(qzero));
    }
  }
  for (size_t i = 0; i < nR; ++i) {
    handlers.rr().fill(0, parameters[i]);
  }
  for (size_t i = 0; i < nR; ++i) {
    for (size_t j = 0; j < nP; ++j)
      handlers.rp().axpy(lin_trans_R(i, j), pparams[j], parameters[i]);
    for (size_t j = 0; j < nQ; ++j)
      handlers.rq().axpy(lin_trans_R(i, nP + j), qparams[j], parameters[i]);
    for (size_t j = 0; j < nD; ++j)
      handlers.rq().axpy(lin_trans_R(i, nP + nQ + j), dparams[j], parameters[i]);
  }
  for (size_t i = 0; i < nDnew; ++i) {
    for (size_t j = 0; j < nP; ++j) {
      handlers.qp().axpy(lin_trans_D(i, j), pparams[j], dparams_new[i]);
      handlers.qp().axpy(lin_trans_D(i, j), pactions[j], dactions_new[i]);
    }
    for (size_t j = 0; j < nQ; ++j) {
      handlers.qq().axpy(lin_trans_D(i, nP + j), qparams[j], dparams_new[i]);
      handlers.qq().axpy(lin_trans_D(i, nP + j), qactions[j], dactions_new[i]);
    }
    for (size_t j = 0; j < nD; ++j) {
      handlers.qq().axpy(lin_trans_D(i, nP + nQ + j), dparams[j], dparams_new[i]);
      handlers.qq().axpy(lin_trans_D(i, nP + nQ + j), dactions[j], dactions_new[i]);
    }
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
  subspace::Matrix<value_type> projected_solutions, overlap;
  std::tie(projected_solutions, overlap) =
      construct_projected_solutions(solutions, xspace.data.at(subspace::EqnData::S), xspace.dimensions());
  auto orth_proj_sol = subspace::Matrix<value_type>{};
  auto norm = std::vector<value_type_abs>{};
  std::tie(orth_proj_sol, norm) = construct_orthogonalised_proj_solutions_in_D(projected_solutions, overlap);
  const auto nC = solutions.rows();
  const auto nD = xspace.dimensions().nD;
  assert(nC >= nD);
  const auto nR = std::min(parameters.size(), xspace.dimensions().nD);
  auto full_solutions = select_full_solutions(orth_proj_sol, norm, nR, norm_thresh_solutions);
  solver.solution_params(full_solutions, parameters);
  const auto nRfull = full_solutions.size();
  const auto params_full = cwrap<R>(begin(parameters), begin(parameters) + nRfull);
  const auto oD = xspace.dimensions().oD;
  auto overlap_DD = subspace::Matrix<value_type>({nD, nD});
  overlap_DD.slice() = xspace.data.at(subspace::EqnData::S).slice({oD, oD}, {oD + nD, oD + nD});
  // FIXME R has components of P+Q subspaces which have to be projected out, else their contribution will be brought
  // into new D
  auto overlap_RD = construct_overlap_RD(params_full, xspace.cparamsd(), overlap_DD, handlers.rr(), handlers.rq());
  auto lin_trans_remaining_D_space = propose_remaining_D_space(overlap_RD, nRfull, norm_thresh_null);
  const auto nRremaining = std::min(nR - nRfull, lin_trans_remaining_D_space.rows());
  auto params_remaining = wrap<R>(begin(parameters) + nRfull, end(parameters));
  construct_remaining_R(params_remaining, lin_trans_remaining_D_space, params_full, xspace.cparamsd(), handlers.rr(),
                        handlers.rq());
  const auto nDnew = lin_trans_remaining_D_space.rows() - nRremaining;
  auto lin_trans_new_D = subspace::Matrix<value_type>({nDnew, nRfull + nD});
  lin_trans_new_D.slice() =
      lin_trans_remaining_D_space.slice({nRfull + nRremaining, 0}, {nRfull + nRremaining + nDnew, nRfull + nD});
  std::vector<Q> dparams, dactions;
  std::tie(dparams, dactions) = construct_new_D(lin_trans_new_D, params_full, xspace.cparamsd(), xspace.cactionsd(),
                                                handlers.qq(), handlers.qr());
  // FIXME might need to normalise the new D parameters and corresponding linear transformation
  auto lin_trans_R_component = subspace::Matrix<value_type>({nDnew, nRfull + nRremaining});
  lin_trans_R_component.slice({0, 0}, {nDnew, nRfull}) = lin_trans_new_D.slice({0, 0}, {nDnew, nRfull});
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
