#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_PROPOSE_RSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_PROPOSE_RSPACE_H
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/subspace/QSpace.h>
#include <molpro/linalg/itsolv/subspace/XSpaceI.h>
#include <molpro/linalg/itsolv/subspace/gram_schmidt.h>
#include <molpro/linalg/itsolv/subspace/util.h>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace detail {

template <class R>
void normalise(VecRef<R>& params, array::ArrayHandler<R, R>& handler, Logger& logger, double thresh = 1.0e-14) {
  for (auto& p : params) {
    auto dot = handler.dot(p, p);
    dot = std::sqrt(std::abs(dot));
    if (dot > thresh) {
      handler.scal(1. / dot, p);
    } else {
      logger.msg("parameter's length is too small for normalisation, dot = " + Logger::scientific(dot), Logger::Warn);
    }
  }
}

//! Proposes an orthonormal set of vectors, removing any params that are linearly dependent
//! @returns linear transformation in orthogonal set and norm
template <class R, typename value_type, typename value_type_abs>
auto propose_orthonormal_set(VecRef<R> params, const double norm_thresh, array::ArrayHandler<R, R>& handler,
                             Logger& logger) {
  logger.msg("propose_orthonormal_set()", Logger::Trace);
  auto lin_trans = subspace::Matrix<value_type>{};
  auto norm = std::vector<value_type_abs>{};
  auto ov = subspace::util::overlap(cwrap(params), handler);
  for (auto done = false; !done;) {
    norm = subspace::util::gram_schmidt(ov, lin_trans);
    auto it_small = std::find_if(begin(norm), end(norm), [norm_thresh](const auto& el) { return el < norm_thresh; });
    done = (it_small == end(norm));
    if (!done) {
      auto i = std::distance(begin(norm), it_small);
      auto it_remove = std::next(begin(params), i);
      params.erase(it_remove);
      ov.remove_row_col(i, i);
      logger.msg("removed parameter index = " + std::to_string(i), Logger::Info);
    }
  }
  lin_trans = subspace::util::construct_lin_trans_in_orthogonal_set(ov, lin_trans, norm);
  return std::tuple<decltype(params), decltype(lin_trans), decltype(norm)>{params, lin_trans, norm};
}
/*!
 * @brief Uses overlap matrix to construct an orthogonal set of R params, and select Q parameters for removal
 *
 *  Linearly dependent Q params are removed based on maximum overlap with proposed R param.
 *
 *  When Q space is larger in size than specified limit, than parameters with smallest maximum contribution to all
 *  solutions are removed.
 *
 * @param overlap overlap matrix of current subspace + initial R parameters
 * @param qspace qspace container
 * @param oQ offset to the Q parameters in the full subspace
 * @param nW number of working parameters
 * @param res_norm_thresh norm threshold for Gram-Schmidt procedure
 * @return index of q parameters to be removed, linear transformation matrix for constructing R params, and their norm
 */
template <typename value_type, typename value_type_abs>
auto calculate_transformation_to_orthogonal_rspace(subspace::Matrix<value_type> overlap,
                                                   const subspace::Matrix<value_type>& solutions,
                                                   const subspace::xspace::Dimensions& dims, Logger& logger,
                                                   value_type_abs res_norm_thresh, unsigned int max_size_qspace) {
  assert(solutions.rows() != 0);
  logger.msg("calculate_transformation_to_orthogonal_rspace()", Logger::Trace);
  auto norm = std::vector<value_type_abs>{};
  auto lin_trans = subspace::Matrix<value_type>{};
  auto qindices_to_remove = std::vector<unsigned int>{};
  auto qindices = std::vector<unsigned int>(dims.nQ);
  std::iota(begin(qindices), end(qindices), 0);
  auto remove_qspace = [&](size_t oQ, size_t iQ) {
    auto iQ_glob = qindices.at(iQ);
    qindices_to_remove.emplace_back(iQ_glob);
    qindices.erase(begin(qindices) + iQ);
    overlap.remove_row_col(oQ + iQ, oQ + iQ);
    logger.msg("removing q parameter = " + std::to_string(iQ_glob), Logger::Info);
  };
  for (bool done = false; !done;) {
    const auto oQ = dims.oQ;
    const auto nQ = qindices.size();
    const auto oN = oQ + nQ;
    norm = subspace::util::gram_schmidt(overlap, lin_trans, res_norm_thresh);
    auto it = std::find_if(std::next(begin(norm), oN), end(norm),
                           [res_norm_thresh](auto el) { return el < res_norm_thresh; });
    auto qspace_is_empty = nQ == 0;
    auto found_singularity = (it == end(norm) && !qspace_is_empty);
    auto qspace_over_limit = nQ > max_size_qspace;
    done = !(found_singularity || qspace_over_limit);
    if (found_singularity) {
      auto i = std::distance(begin(norm), it);
      logger.msg("parameter index i = " + std::to_string(i) + " norm = " + std::to_string(*it), Logger::Info);
      auto normalised_overlap = std::vector<value_type>{};
      for (size_t j = 0; j < nQ; ++j)
        normalised_overlap.emplace_back(std::abs(overlap(i, oQ + j)) / std::sqrt(std::abs(overlap(oQ + j, oQ + j))));
      auto it_max = std::max_element(begin(normalised_overlap), end(normalised_overlap));
      auto iq_erase = std::distance(begin(normalised_overlap), it_max);
      remove_qspace(oQ, iq_erase);
    } else if (qspace_over_limit) {
      auto max_contrib_to_solution = std::vector<value_type_abs>{};
      for (auto i : qindices) {
        const auto nSol = solutions.rows();
        auto contrib = std::vector<value_type_abs>(nSol);
        for (size_t j = 0; j < nSol; ++j) {
          contrib[j] = std::abs(solutions(j, oQ + i));
        }
        max_contrib_to_solution.emplace_back(*std::max_element(begin(contrib), end(contrib)));
      }
      auto it_min = std::min_element(begin(max_contrib_to_solution), end(max_contrib_to_solution));
      auto i = std::distance(begin(max_contrib_to_solution), it_min);
      remove_qspace(oQ, i);
    }
  }
  return std::tuple<decltype(qindices_to_remove), decltype(lin_trans), decltype(norm)>{qindices_to_remove, lin_trans,
                                                                                       norm};
}

/*!
 * @brief Projects solution from the full subspace on to Q_{delete} and current D space.
 * @param solutions solution matrix in the full subspace
 * @param dims dimensions for partitioning of subspace
 * @param remove_qspace indices to remove from current Q space and move into Q_{delete}
 * @param overlap overlap matrix of the current subspace
 * @param norm_thresh vectors with norm less than this threshold are considered null
 * @return
 */
template <typename value_type, typename value_type_abs>
auto construct_projected_solution(const subspace::Matrix<value_type>& solutions,
                                  const subspace::xspace::Dimensions& dims,
                                  const std::vector<unsigned int>& remove_qspace,
                                  const subspace::Matrix<value_type>& overlap, value_type_abs norm_thresh) {
  const auto nQd = remove_qspace.size();
  const auto nSol = solutions.rows();
  auto solutions_proj = subspace::Matrix<value_type>({nSol, nQd + dims.nD});
  for (size_t i = 0; i < nSol; ++i) {
    for (size_t j = 0; j < nQd; ++j) {
      solutions_proj(i, j) = solutions(i, remove_qspace[j]);
    }
    for (size_t j = 0; j < dims.nD; ++j) {
      solutions_proj(i, nQd + j) = solutions(i, dims.oD + j);
    }
  }
  auto norm_proj = std::vector<value_type_abs>(nSol);
  for (size_t i = 0; i < nSol; ++i) {
    for (size_t j = 0; j < nQd; ++j) {
      for (size_t k = 0; k < j; ++k) {
        norm_proj[i] += 2. * solutions_proj(i, j) * solutions_proj(i, k) * overlap(remove_qspace[j], remove_qspace[k]);
      }
      norm_proj[i] += std::pow(solutions_proj(i, j), 2) * overlap(remove_qspace[j], remove_qspace[j]);
    }
    for (size_t j = 0; j < nQd; ++j) {
      for (size_t k = 0; k < dims.nD; ++k) {
        norm_proj[i] += 2. * solutions_proj(i, j) * solutions_proj(i, nQd + k) * overlap(remove_qspace[j], dims.oD + k);
      }
    }
    for (size_t j = 0; j < dims.nD; ++j) {
      for (size_t k = 0; k < j; ++k) {
        norm_proj[i] += 2. * solutions_proj(nQd + i, nQd + j) * solutions_proj(nQd + i, nQd + k) *
                        overlap(dims.oD + j, dims.oD + k);
      }
      norm_proj[i] += std::pow(solutions_proj(nQd + i, nQd + j), 2) * overlap(dims.oD + j, dims.oD + j);
    }
    for (size_t j = 0; j < dims.nD; ++j) {
      solutions_proj(i, nQd + j) = solutions(i, dims.oD + j);
    }
  }
  for (auto& x : norm_proj)
    x = std::sqrt(std::abs(x));
  for (size_t i = 0, j = 0; i < nSol; ++i) {
    if (norm_proj[i] > norm_thresh) {
      solutions_proj.row(j).scal(1. / norm_proj[i]);
      ++j;
    } else {
      solutions_proj.remove_row(j);
    }
  }
  return solutions_proj;
}

/*!
 * @brief Constructs overlap matrix of P+Q+R+(projected solutions) subspaces, where Q is without removed parameters
 * @param solutions_proj solutions matrix with Q deleted + current D space as columns
 * @param dims dimensions of the current subspace
 * @param remove_qspace indices of Q parameters to be removed
 * @param overlap overlap of current subspace P+Q+D+R, including all of the current Q space
 * @param nR number of new parameters
 * @return overlap matrix for the subspace P+Q+R+(projected solutions), Q without deleted vectors
 */
template <typename value_type>
auto construct_overlap_with_projected_solutions(const subspace::Matrix<value_type>& solutions_proj,
                                                const subspace::xspace::Dimensions& dims,
                                                const std::vector<unsigned int>& remove_qspace,
                                                const subspace::Matrix<value_type>& overlap, const size_t nR) {
  const auto nDnew = solutions_proj.size();
  const auto nQd = remove_qspace.size();
  const auto nQ = dims.nQ - nQd;
  auto ov = overlap;
  for (size_t i = 0; i < dims.nD; ++i) {
    ov.remove_row_col(dims.oD, dims.oD);
  }
  auto is_Qdelete = [&remove_qspace](size_t i) {
    return std::find(begin(remove_qspace), end(remove_qspace), i) != end(remove_qspace);
  };
  for (size_t i = 0, j = 0; i < dims.nQ; ++i) {
    if (is_Qdelete(i))
      ov.remove_row_col(dims.oQ + j, dims.oQ + j);
    else
      ++j;
  }
  const auto oDnew = dims.nP + nQ + nR;
  ov.resize({oDnew + nDnew, oDnew + nDnew});
  /*
   * append overlap with solutions_proj
   * x_i = \sum_j c_ij v_j
   * <x_i, v_j> = \sum_k c_ik <v_k, v_j>
   * <x_i, x_j> = \sum_kl c_ik c_jl <v_k, v_l>
   */
  auto accumulate_ov_offdiag = [&](size_t i, size_t j, size_t jj) {
    for (size_t k = 0; k < nQd; ++k)
      ov(oDnew + i, j) += solutions_proj(i, k) * overlap(jj, dims.oQ + remove_qspace[k]);
    for (size_t k = 0; k < dims.nD; ++k)
      ov(oDnew + i, j) += solutions_proj(i, nQd + k) * overlap(jj, dims.oD + k);
    ov(j, oDnew + i) = ov(oDnew + i, j);
  };
  for (size_t i = 0; i < nDnew; ++i) {
    for (size_t j = 0; j < dims.nP; ++j)
      accumulate_ov_offdiag(i, j, dims.oP + j);
    for (size_t j = 0, jj = 0; j < dims.nQ; ++j)
      if (!is_Qdelete(j))
        accumulate_ov_offdiag(i, jj++, dims.oQ + j);
    for (size_t j = 0; j < nR; ++j)
      accumulate_ov_offdiag(i, dims.nP + nQ + j, dims.nX + j);
    for (size_t j = 0; j <= i; ++j) {
      for (size_t k = 0; k < nQd; ++k) {
        for (size_t l = 0; l < nQd; ++l)
          ov(oDnew + i, oDnew + j) += solutions_proj(i, k) * solutions_proj(j, l) *
                                      overlap(dims.oQ + remove_qspace[k], dims.oQ + remove_qspace[l]);
        for (size_t l = 0; l < dims.nD; ++l)
          ov(oDnew + i, oDnew + j) +=
              solutions_proj(i, k) * solutions_proj(j, nQd + l) * overlap(dims.oQ + remove_qspace[k], dims.oD + l);
      }
      for (size_t k = 0; k < dims.nD; ++k) {
        for (size_t l = 0; l < nQd; ++l)
          ov(oDnew + i, oDnew + j) +=
              solutions_proj(i, nQd + k) * solutions_proj(j, l) * overlap(dims.oD + k, dims.oQ + remove_qspace[l]);
        for (size_t l = 0; l < dims.nD; ++l)
          ov(oDnew + i, oDnew + j) +=
              solutions_proj(i, nQd + k) * solutions_proj(j, nQd + l) * overlap(dims.oD + k, dims.oD + l);
      }
      ov(oDnew + j, oDnew + i) = ov(oDnew + i, oDnew + j);
    }
  }
  return ov;
}

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
 * @brief Constructs transformation to D space by projecting solutions on to deleted Q (Q_D) and old D space, than
 *  orthogonalising against P+Q+R
 * @param solutions row-wise matrix with solutions in the subspace
 * @param dims dimensions of the current subspace
 * @param remove_qspace indices of Q parameters that are deleted in this iteration
 * @param overlap overlap matrix for current subspace with new R parameters (P + Q + D + R)
 * @return linear transformation to the new D space in terms of P+Q+R+Qdelete+Dold
 */
template <typename value_type, typename value_type_abs>
auto propose_dspace(const subspace::Matrix<value_type>& solutions, const subspace::xspace::Dimensions& dims,
                    const std::vector<unsigned int>& remove_qspace, const subspace::Matrix<value_type>& overlap,
                    const size_t nR, value_type_abs norm_thresh) {
  auto solutions_proj = construct_projected_solution(solutions, dims, remove_qspace, overlap, norm_thresh);
  // construct overlap of projected solutions with the new P+Q+R subspace (excluding Q_D)
  auto ov = construct_overlap_with_projected_solutions(solutions_proj, dims, remove_qspace, overlap, nR);
  // orthogonalise against the subspace
  auto lin_trans = subspace::Matrix<value_type>{};
  auto norm = subspace::util::gram_schmidt(ov, lin_trans);
  const auto nSol = solutions_proj.size();
  const auto nX = norm.size() - nSol;
  remove_null_vectors(lin_trans, norm, nX, norm.size(), norm_thresh);
  const auto nD = norm.size() - nX;
  // Orthogonalised D space is in terms of projected solutions
  // I need to use the old D space instead
  const auto nY = remove_qspace.size() + dims.nD;
  auto lin_trans_Dold = subspace::Matrix<value_type>({nD, nX + nY});
  lin_trans_Dold.slice({0, 0}, {nD, nX}) = lin_trans.slice({nX, 0}, {nX + nD, nX});
  /*
   * x_i = \sum_j T_ij v_j
   * v_i = \sum_j C_ij u_j
   * x_i = \sum_j \sum_k T_ij C_jk u_k
   * x_i = \sum_k (\sum_j T_ij C_jk) u_k
   * x_i = \sum_k D_ik u_k
   * D_ij = \sum_k T_ik C_kj
   */
  for (size_t i = 0; i < nD; ++i)
    for (size_t j = 0; j < nY; ++j)
      for (size_t k = 0; k < nSol; ++k)
        lin_trans_Dold(i, j) += lin_trans(i, k) * solutions_proj(k, j);
  norm.erase(begin(norm), begin(norm) + nX);
  return std::tuple<decltype(lin_trans), decltype(norm)>{lin_trans, norm};
}

/*!
 * @brief Applies linear transformation to construct the D space
 * @param xspace subspace container. New D space is stored in xspace directly
 * @param lin_trans new D space vectors in the subspace (P+Q+R+Qdelete+Dold)
 * @param norm estimated norm of transformed vectors
 * @param rparams R space parameters
 * @param handlers array handlers
 */
template <class R, class Q, class P, typename value_type, typename value_type_abs>
void construct_orthonormal_Dparams(subspace::XSpaceI<R, Q, P>& xspace, const subspace::Matrix<value_type>& lin_trans,
                                   const std::vector<value_type_abs>& norm, const CVecRef<Q>& rparams,
                                   ArrayHandlers<R, Q, P>& handlers) {}

/*!
 * @brief Construct an orthonormal set from params and a linear transformation matrix
 *
 * @param params vectors to orthonormalise
 * @param lin_trans Gram Schmidt linear transformation in an orthogonal set, see construct_lin_trans_in_orthogonal_set()
 * @param norm estimated norm of orthogonalised vectors
 */
template <class R, typename value_type, typename value_type_abs>
void construct_orthonormal_set(VecRef<R>& params, const subspace::Matrix<value_type>& lin_trans,
                               const std::vector<value_type_abs>& norm, array::ArrayHandler<R, R>& handler) {
  for (size_t i = 0; i < params.size(); ++i) {
    for (size_t j = 0; j < i; ++j) {
      handler.axpy(lin_trans(i, j), params.at(j), params.at(i));
    }
  }
  for (size_t i = 0; i < params.size(); ++i) {
    handler.scal(1. / norm.at(i), params.at(i));
  }
}
/*!
 * @brief Appends row and column for overlap with params
 *
 * Overlap with previous solutions is not included.
 *
 * @param ov overlap matrix
 * @param params new parameters
 * @param pparams P space parameters
 * @param qparams Q space parameters
 * @param handlers array handlers
 * @param logger logger
 */
template <class R, class Q, class P, typename value_type>
auto append_overlap_with_r(const subspace::Matrix<value_type>& overlap, const CVecRef<R>& params,
                           const CVecRef<P>& pparams, const CVecRef<Q>& qparams, ArrayHandlers<R, Q, P>& handlers,
                           Logger& logger) {
  auto nP = pparams.size();
  auto nQ = qparams.size();
  auto nN = params.size();
  auto nX = nP + nQ + nN;
  auto oP = 0;
  auto oQ = oP + nP;
  auto oN = oQ + nQ;
  auto ov = overlap;
  ov.resize({nX, nX}); // solutions are last
  ov.slice({oN, oN}, {oN + nN, oN + nN}) = subspace::util::overlap(params, handlers.rr());
  ov.slice({oN, oP}, {oN + nN, oP + nP}) = subspace::util::overlap(params, pparams, handlers.rp());
  ov.slice({oN, oQ}, {oN + nN, oQ + nQ}) = subspace::util::overlap(params, qparams, handlers.rq());
  auto copy_upper_to_lower = [&ov, oN, nN](size_t oX, size_t nX) {
    for (size_t i = 0; i < nX; ++i)
      for (size_t j = 0; j < nN; ++j)
        ov(oX + i, oN + j) = ov(oN + j, oX + i);
  };
  copy_upper_to_lower(oP, nP);
  copy_upper_to_lower(oQ, nQ);
  return ov;
}
/*!
 * @brief Constructs R parameters approximately orthonormal in the full subspace
 * @param params output new R parameters
 * @param residuals input non-orthogonal parameters
 * @param lin_trans linear transformation in the full subspace
 * @param norm approximate norm after linear transformation (full subspace)
 * @param pparams P space parameters
 * @param qparams Q space parameters
 * @param handlers
 */
template <class R, class Q, class P, typename value_type, typename value_type_abs>
void construct_orthonormal_Rparams(VecRef<R>& params, VecRef<R>& residuals,
                                   const subspace::Matrix<value_type>& lin_trans,
                                   const std::vector<value_type_abs>& norm, const CVecRef<P>& pparams,
                                   const CVecRef<Q>& qparams, ArrayHandlers<R, Q, P>& handlers) {
  assert(params.size() == residuals.size());
  auto nP = pparams.size();
  auto nQ = qparams.size();
  auto nN = params.size();
  auto oP = 0;
  auto oQ = oP + nP;
  auto oN = oQ + nQ;
  for (size_t i = 0; i < nN; ++i) {
    handlers.rr().copy(params.at(i), residuals.at(i));
  }
  for (size_t i = 0; i < nN; ++i) {
    for (size_t j = 0; j < nP; ++j) {
      handlers.rp().axpy(lin_trans(oN + i, oP + j), pparams.at(j), params.at(i));
    }
    for (size_t j = 0; j < nQ; ++j) {
      handlers.rq().axpy(lin_trans(oN + i, oQ + j), qparams.at(j), params.at(i));
    }
    for (size_t j = 0; j < i; ++j) {
      handlers.rr().axpy(lin_trans(oN + i, oN + j), residuals.at(j), params.at(i));
    }
  }
  for (size_t i = 0; i < nN; ++i) {
    handlers.rr().scal(1. / norm.at(oN + i), params.at(i));
  }
}

//! Returns new working set based on parameters included in wparams
template <class R>
auto get_new_working_set(const std::vector<unsigned int>& working_set, const std::vector<R>& params,
                         const VecRef<R>& wparams) {
  auto new_indices = find_ref(wparams, begin(params), end(params));
  auto new_working_set = std::vector<unsigned int>{};
  for (auto i : new_indices) {
    new_working_set.emplace_back(working_set.at(i));
  }
  return new_working_set;
}

/*!
 * \brief Proposes new parameters for the subspace from the preconditioned residuals.
 *
 * Outline
 * -------
 * Basic procedure:
 *  - Gram-schmidt orthonormalise residuals amongst themselves
 *  - Gram-schmidt orthogonalise residuals against the old Q space and current solutions
 *  - Ensure that the resultant Q space is not linearly dependent by removing Q parameters with large overlap
 *  - Ensure that the size of Q space is within limit, by removing Q parameters with smallest contributions to current
 *    solutions
 *  - Reconstruct C space
 *
 * Various possibilities:
 *  1. Residuals are linearly dependent among themselves, worst case scenario there could be duplicates.
 *  2. Residuals are linearly dependent with the old Q space, orthonormalisation against Q would result in
 *     almost null vectors.
 *
 * Case 1 is handled at the start by normalising residuals and orthogonalising them among themselves. If it results in
 * vectors with norm less then **threshold** than they are discarded and their action does not need to be evaluated.
 *
 * Case 2 is handled during Gram-Schmidt procedure. Residuals are orthogonalised against the old Q space, if one of
 * them has a small norm than an old q vector with largest overlap is deleted.
 *
 * @param parameters output new parameters for the subspace.
 * @param residual preconditioned residuals.
 * @return number of significant parameters to calculate the action for
 */
template <class R, class Q, class P, typename value_type, typename value_type_abs>
auto propose_rspace(LinearEigensystem<R, Q, P>& solver, std::vector<R>& parameters, std::vector<R>& residuals,
                    subspace::XSpaceI<R, Q, P>& xspace, const subspace::Matrix<value_type>& solutions,
                    ArrayHandlers<R, Q, P>& handlers, Logger& logger, value_type_abs res_norm_thresh,
                    unsigned int max_size_qspace) {
  logger.msg("itsolv::detail::propose_rspace", Logger::Trace);
  auto wresidual = wrap<R>(residuals.begin(), residuals.begin() + solver.working_set().size());
  normalise(wresidual, handlers.rr(), logger);
  auto lin_trans = subspace::Matrix<value_type>{};
  auto norm = std::vector<value_type_abs>{};
  std::tie(wresidual, lin_trans, norm) =
      propose_orthonormal_set<R, value_type, value_type_abs>(wresidual, res_norm_thresh, handlers.rr(), logger);
  // FIXME propose and construct in one go
  construct_orthonormal_set(wresidual, lin_trans, norm, handlers.rr());
  // propose working space by orthogonalising against P+Q
  auto ov = append_overlap_with_r(xspace.data.at(subspace::EqnData::S), cwrap(wresidual), xspace.cparamsp(),
                                  xspace.cparamsq(), handlers, logger);
  auto remove_qspace = std::vector<unsigned int>{};
  std::tie(remove_qspace, lin_trans, norm) = calculate_transformation_to_orthogonal_rspace(
      ov, solutions, xspace.dimensions(), logger, res_norm_thresh, max_size_qspace);
  if (logger.data_dump) {
    logger.msg("full overlap = " + subspace::as_string(ov), Logger::Info);
    logger.msg("linear transformation = " + subspace::as_string(lin_trans), Logger::Info);
    logger.msg("norm = ", norm.begin(), norm.end(), Logger::Info);
    logger.msg("remove Q space indices = ", remove_qspace.begin(), remove_qspace.end(), Logger::Info);
  }
  auto new_working_set = get_new_working_set(solver.working_set(), residuals, wresidual);
  auto wparams = wrap<R>(parameters.begin(), parameters.begin() + wresidual.size());
  auto qparams_new = remove_elements(xspace.cparamsq(), remove_qspace);
  construct_orthonormal_Rparams(wparams, wresidual, lin_trans, norm, xspace.cparamsp(), qparams_new, handlers);
  normalise(wparams, handlers.rr(), logger);
  auto params_qd = xspace.cparamsq();
  auto paramsd = xspace.cparamsd();
  std::copy(begin(paramsd), end(paramsd), std::back_inserter(params_qd));
  ov = append_overlap_with_r(xspace.data.at(subspace::EqnData::S), cwrap(wparams), xspace.cparamsp(), params_qd,
                             handlers, logger);
  std::tie(lin_trans, norm) =
      propose_dspace(solutions, xspace.dimensions(), remove_qspace, ov, wparams.size(), res_norm_thresh);
  construct_orthonormal_Dparams(xspace, lin_trans, norm, remove_qspace, cwrap(wparams), handlers);
  // finally remove Q parameters
  return new_working_set;
}
} // namespace detail
} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_PROPOSE_RSPACE_H
