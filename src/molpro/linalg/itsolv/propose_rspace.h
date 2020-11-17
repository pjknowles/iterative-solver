#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_PROPOSE_RSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_PROPOSE_RSPACE_H
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/helper.h>
#include <molpro/linalg/itsolv/subspace/Dimensions.h>
#include <molpro/linalg/itsolv/subspace/QSpace.h>
#include <molpro/linalg/itsolv/subspace/XSpaceI.h>
#include <molpro/linalg/itsolv/subspace/gram_schmidt.h>
#include <molpro/linalg/itsolv/subspace/util.h>
#include <molpro/linalg/itsolv/util.h>
#include <molpro/linalg/itsolv/wrap_util.h>

namespace molpro::linalg::itsolv::detail {

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
  return std::make_tuple(params, lin_trans, norm);
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
 * @param svd_thresh theshold on singular values to mark that the subspace is not stable
 * @param norm_thresh norm threshold for Gram-Schmidt procedure
 * @return index of q parameters to be removed, linear transformation matrix for constructing R params, and their norm
 */
template <typename value_type, typename value_type_abs>
auto calculate_transformation_to_orthogonal_rspace(subspace::Matrix<value_type> overlap,
                                                   const subspace::Matrix<value_type>& solutions,
                                                   const subspace::Dimensions& dims, Logger& logger,
                                                   value_type_abs svd_thresh, value_type_abs norm_thresh,
                                                   int max_size_qspace) {
  assert(solutions.rows() != 0);
  logger.msg("calculate_transformation_to_orthogonal_rspace()", Logger::Trace);
  const auto oQ = dims.oQ;
  auto norm = std::vector<value_type_abs>{};
  auto lin_trans = subspace::Matrix<value_type>{};
  auto qindices_to_remove = std::vector<int>{};
  auto qindices = std::vector<int>(dims.nQ);
  std::iota(begin(qindices), end(qindices), 0);
  auto remove_qspace = [&](size_t iQ) {
    auto iQ_glob = qindices.at(iQ);
    qindices_to_remove.emplace_back(iQ_glob);
    qindices.erase(begin(qindices) + iQ);
    overlap.remove_row_col(oQ + iQ, oQ + iQ);
    logger.msg("removing q parameter = " + std::to_string(iQ_glob), Logger::Info);
  };
  for (bool found_singularity = true; found_singularity;) {
    const auto nQ = qindices.size();
    auto svd_sys = svd_system(overlap.rows(), overlap.cols(), array::Span(&overlap(0, 0), overlap.size()), svd_thresh);
    auto qspace_is_empty = nQ == 0;
    found_singularity = (!svd_sys.empty() && !qspace_is_empty);
    if (found_singularity) {
      auto sing_val_from_q = std::vector<value_type_abs>(nQ);
      for (size_t i = 0; i < nQ; ++i)
        sing_val_from_q[i] = std::abs(svd_sys.front().v.at(oQ + i));
      auto it_max = std::max_element(begin(sing_val_from_q), end(sing_val_from_q));
      auto iq_erase = std::distance(begin(sing_val_from_q), it_max);
      {
        std::stringstream sv;
        sv << std::setprecision(3) << svd_sys.front().value;
        logger.msg("found singularity due to Q parameter index = " + std::to_string(iq_erase) + ", singular value = " +
                       sv.str() + ", SVD right vector contribution = " + std::to_string(*it_max),
                   Logger::Debug);
      }
      remove_qspace(iq_erase);
    }
  }
  while (qindices.size() > max_size_qspace) {
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
    remove_qspace(i);
  }
  norm = subspace::util::gram_schmidt(overlap, lin_trans, norm_thresh);
  return std::make_tuple(qindices_to_remove, lin_trans, norm);
}
namespace dspace {

/*!
 * @brief Projects solution from the full subspace on to Q_{delete} and current D space.
 * @param solutions solution matrix in the full subspace
 * @param dims dimensions for partitioning of subspace
 * @param remove_qspace indices to remove from current Q space and move into Q_{delete}
 * @param norm_thresh vectors with norm less than this threshold are considered null
 * @return projected solutions, and their overlap matrix
 */
template <typename value_type>
auto construct_projected_solution(const subspace::Matrix<value_type>& solutions, const subspace::Dimensions& dims,
                                  const std::vector<int>& remove_qspace, Logger& logger) {
  logger.msg("construct_projected_solution()", Logger::Trace);
  const auto nQd = remove_qspace.size();
  const auto nSol = solutions.rows();
  auto solutions_proj = subspace::Matrix<value_type>({nSol, nQd + dims.nD});
  for (size_t i = 0; i < nSol; ++i) {
    for (size_t j = 0; j < nQd; ++j) {
      solutions_proj(i, j) = solutions(i, dims.oQ + remove_qspace[j]);
    }
    for (size_t j = 0; j < dims.nD; ++j) {
      solutions_proj(i, nQd + j) = solutions(i, dims.oD + j);
    }
  }
  logger.msg("nSol, nQd, nD " + std::to_string(nSol) + ", " + std::to_string(nQd) + ", " + std::to_string(dims.nD),
             Logger::Debug);
  return solutions_proj;
}

/*!
 * @brief Constructs overlap matrix for projected solutions
 *
 * x_i = \sum_j C_ij u_j
 * <x_i, x_i> = \sum_j \sum_k C_ij C_ik <u_j, u_k>
 * <x_i, x_j> = \sum_k \sum_l C_ik C_jl <u_j, u_l>
 *
 * @param solutions_proj solutions matrix projected onto Qd+D space
 * @param dims dimensions for partitioning of subspace
 * @param remove_qspace indices to remove from current Q space and move into Q_{delete}
 * @param overlap overlap matrix of the current subspace
 * @param logger logger
 * @return overlap matrix
 */
template <typename value_type>
auto construct_projected_solutions_overlap(const subspace::Matrix<value_type>& solutions_proj,
                                           const subspace::Matrix<value_type>& overlap,
                                           const subspace::Dimensions& dims, const std::vector<int>& remove_qspace,
                                           Logger& logger) {
  logger.msg("construct_projected_solution_overlap()", Logger::Trace);
  const auto nSol = solutions_proj.rows();
  const auto nQd = remove_qspace.size();
  auto overlap_proj = subspace::Matrix<value_type>({nSol, nSol});
  for (size_t i = 0; i < nSol; ++i) {
    for (size_t ii = 0; ii <= i; ++ii) {
      for (size_t j = 0; j < nQd; ++j) {
        for (size_t k = 0; k < nQd; ++k) {
          overlap_proj(i, ii) += solutions_proj(i, j) * solutions_proj(ii, k) *
                                 overlap(dims.oQ + remove_qspace[j], dims.oQ + remove_qspace[k]);
        }
        for (size_t k = 0; k < dims.nD; ++k) {
          overlap_proj(i, ii) +=
              solutions_proj(i, j) * solutions_proj(ii, nQd + k) * overlap(dims.oQ + remove_qspace[j], dims.oD + k);
        }
      }
      for (size_t j = 0; j < dims.nD; ++j) {
        for (size_t k = 0; k < dims.nD; ++k) {
          overlap_proj(i, ii) +=
              solutions_proj(i, nQd + j) * solutions_proj(ii, nQd + k) * overlap(dims.oD + j, dims.oD + k);
        }
        for (size_t k = 0; k < nQd; ++k) {
          overlap_proj(i, ii) += solutions_proj(i, nQd + j) * solutions_proj(ii, k) * overlap(dims.oD + j, dims.oQ + k);
        }
      }
      overlap_proj(ii, i) = overlap_proj(i, ii);
    }
  }
  return overlap_proj;
}

/*!
 * @brief Removes parameters with norm less than threshold and normalises the rest
 * @param parameters row major matrix of parameters
 * @param overlap overlap matrix of parameters
 * @param norm_thresh norm threhsold
 * @param logger logger
 */
template <typename value_type, typename value_type_abs>
void remove_null_norm_and_normalise(subspace::Matrix<value_type>& parameters, subspace::Matrix<value_type>& overlap,
                                    const value_type_abs norm_thresh, Logger& logger) {
  logger.msg("remove_null_norm_and_normalise()", Logger::Trace);
  const auto nSol = parameters.rows();
  auto norm_proj = std::vector<value_type_abs>(nSol, 0.);
  for (size_t i = 0; i < nSol; ++i)
    norm_proj[i] = std::sqrt(std::abs(overlap(i, i)));
  for (size_t i = 0, j = 0; i < nSol; ++i) {
    if (norm_proj[i] > norm_thresh) {
      parameters.row(j).scal(1. / norm_proj[i]);
      overlap.col(j).scal(1. / norm_proj[i]);
      overlap.row(j).scal(1. / norm_proj[i]);
      ++j;
    } else {
      parameters.remove_row(j);
      overlap.remove_row_col(j, j);
    }
  }
  if (logger.data_dump) {
    logger.msg("norm  = ", std::begin(norm_proj), std::end(norm_proj), Logger::Info);
    logger.msg("parameters after normalisation = " + as_string(parameters), Logger::Info);
    logger.msg("overlap of parameters = " + as_string(overlap), Logger::Info);
  }
}

/*!
 * @brief Normalises projected solutions and removes the null space if nSol > nQd+nD
 * @param solutions_proj solutions porjected onto Qd+D space
 * @param solutions solutions in the current subspace
 * @param overlap_proj overlap matrix of projected solutions
 * @param norm_thresh
 * @param logger
 */
template <typename value_type, typename value_type_abs>
auto remove_null_projected_solutions(subspace::Matrix<value_type>& solutions_proj,
                                     subspace::Matrix<value_type>& overlap_proj,
                                     const subspace::Matrix<value_type>& solutions, const value_type_abs norm_thresh,
                                     const value_type_abs svd_thresh, Logger& logger) {
  auto svd_vecs =
      svd_system(overlap_proj.rows(), overlap_proj.cols(), array::Span(&overlap_proj(0, 0), overlap_proj.size()),
                 std::numeric_limits<value_type_abs>::max());
  std::remove_if(std::begin(svd_vecs), std::end(svd_vecs),
                 [&svd_thresh](const auto& el) { return el.value < svd_thresh; });
  size_t n_removed = 0;
  while (solutions_proj.rows() > solutions_proj.cols()) {
    auto lin_trans = subspace::Matrix<value_type>{};
    auto norm = subspace::util::gram_schmidt(overlap_proj, lin_trans, norm_thresh);
    auto it = std::min_element(std::begin(norm), std::end(norm));
    auto i = std::distance(std::begin(norm), it);
    solutions_proj.remove_row(i);
    solutions_D.remove_row(i);
    overlap_proj.remove_row_col(i, i);
    logger.msg("removing null projected solution, i = " + std::to_string(i) + ", norm = " + std::to_string(*it),
               Logger::Debug);
    ++n_removed;
  }
  logger.msg("removed n = " + std::to_string(n_removed) + " projected solutions", Logger::Debug);
  logger.msg("remaining number of projected solutions = " + std::to_string(solutions_proj.rows()), Logger::Debug);
  svd_vecs =
      svd_system(overlap_proj.rows(), overlap_proj.cols(), array::Span(&overlap_proj(0, 0), overlap_proj.size()), 1);
  return solutions_D;
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
                                                const subspace::Dimensions& dims, const std::vector<int>& remove_qspace,
                                                const subspace::Matrix<value_type>& overlap, const size_t nR) {
  const auto nDnew = solutions_proj.rows();
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
        accumulate_ov_offdiag(i, dims.nP + jj++, dims.oQ + j);
    for (size_t j = 0; j < nR; ++j)
      accumulate_ov_offdiag(i, dims.nP + nQ + j, dims.nX + j);
  }
  for (size_t i = 0; i < nDnew; ++i) {
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
 * @brief Constructs transformation to D space by projecting solutions on to deleted Q (Q_D) and old D space, than
 *  orthogonalising against P+Q+R
 * @param solutions row-wise matrix with solutions in the subspace
 * @param dims dimensions of the current subspace
 * @param remove_qspace indices of Q parameters that are deleted in this iteration
 * @param overlap overlap matrix for current subspace with new R parameters (P + Q + D + R)
 * @return linear transformation to the new D space in terms of P+Q+R+Qdelete+Dold, and corresponding solutions
 */
template <typename value_type, typename value_type_abs>
auto propose_dspace(const subspace::Matrix<value_type>& solutions, const subspace::Dimensions& dims,
                    std::vector<int>& remove_qspace, const subspace::Matrix<value_type>& overlap, const size_t nR,
                    value_type_abs norm_thresh, Logger& logger) {
  logger.msg("propose_dspace()", Logger::Trace);
  auto solutions_proj = construct_projected_solution(solutions, dims, remove_qspace, logger);
  auto overlap_proj = construct_projected_solutions_overlap(solutions_proj, overlap, dims, remove_qspace, logger);
  remove_null_norm_and_normalise(solutions_proj, overlap_proj, norm_thresh, logger);
  // perform SVD and remove null space by transforming to the basis of significant right-singular vectors
  auto solutions_D =
      remove_null_projected_solutions(solutions_proj, overlap_proj, solutions, norm_thresh, norm_thresh, logger);
  // There is no mapping from D parameters to individual solutions
  auto ov = construct_overlap_with_projected_solutions(solutions_proj, dims, remove_qspace, overlap, nR);
  // Do SVD of the full subspace.
  // If there are null parameters, remove Q with smallest contribution to the solution
  // repeat until the subspace is well conditioned
  // If there are no Q parameters left, remove D parameters with largest contributions to the null space
  auto svd_vecs = svd_system(ov.rows(), ov.cols(), array::Span(&ov(0, 0), ov.size()), 1);
  auto lin_trans = subspace::Matrix<value_type>{};
  auto norm = subspace::util::gram_schmidt(ov, lin_trans);
  if (logger.data_dump) {
    logger.msg("overlap matrix of P+Q+R+projected solution = " + as_string(ov), Logger::Info);
    logger.msg("lin_trans = " + as_string(lin_trans), Logger::Info);
    logger.msg("norm = ", begin(norm), end(norm), Logger::Info);
  }
  const auto nSol = solutions_proj.rows();
  const auto nX = norm.size() - nSol;
  {
    auto norm_D = std::vector<value_type_abs>(norm.begin() + nX, norm.end());
    util::remove_null_vectors(lin_trans, norm, nX, norm.size(), norm_thresh);
    util::remove_null_vectors(solutions_D, norm_D, 0, norm_D.size(), norm_thresh);
  }
  auto nD = norm.size() - nX;
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
        lin_trans_Dold(i, nX + j) += lin_trans(nX + i, nX + k) * solutions_proj(k, j);
  for (size_t i = 0; i < nD; ++i)
    lin_trans_Dold.row(i).scal(1. / norm[nX + i]);
  const auto nDmax = std::min(nSol, remove_qspace.size() + dims.nD);
  while (nD > nDmax) {
    auto it = std::min_element(std::begin(norm) + nX, std::end(norm));
    auto i = std::distance(std::begin(norm) + nX, it);
    logger.msg("nD > nDmax, {nD, nDmax} = " + std::to_string(nD) + ", " + std::to_string(nDmax), Logger::Debug);
    logger.msg("erase i, norm = " + std::to_string(i) + ", " + std::to_string(*it), Logger::Debug);
    norm.erase(it);
    lin_trans_Dold.remove_row(i);
    solutions_D.remove_row(i);
    --nD;
  }
  logger.msg("nD = " + std::to_string(nD), Logger::Debug);
  if (logger.data_dump) {
    logger.msg("lin_trans using old D = " + as_string(lin_trans_Dold), Logger::Info);
    logger.msg("solutions D = " + as_string(solutions_D), Logger::Info);
  }
  return std::make_tuple(lin_trans_Dold, solutions_D);
}
} // namespace dspace

/*!
 * @brief Applies linear transformation to construct the D space parameters and corresponding action (without R space
 * action, since it is not yet known)
 *
 * Subsystem Hamiltonian in D space and D space action cannot be constructed until R space action is known.
 * Thus, it becomes simpler to calculate the subspace equation data using overlap of full vectors instead of subspace
 * arithmetic.
 *
 * @param xspace subspace container. New D space is stored in xspace directly
 * @param lin_trans new D space vectors in the subspace (P+Q+R+Qdelete+Dold)
 * @param rparams R space parameters
 * @param handlers array handlers
 * @returns D space parameters, actions and R space component of the linear transformation matrix
 */
template <class R, class Q, class P, typename value_type>
auto construct_orthonormal_Dparams(subspace::XSpaceI<R, Q, P>& xspace, const subspace::Matrix<value_type>& lin_trans,
                                   const std::vector<int>& q_indices_remove, const CVecRef<R>& rparams,
                                   ArrayHandlers<R, Q, P>& handlers, Logger& logger) {
  const auto nD = lin_trans.rows();
  const auto nQdelete = q_indices_remove.size();
  const auto qparams = xspace.cparamsq();
  const auto pparams = xspace.cparamsp();
  const auto dparams_old = xspace.cparamsd();
  const auto qactions = xspace.cactionsq();
  const auto pactions = xspace.cactionsp();
  const auto dactions_old = xspace.cactionsd();
  const auto dims = xspace.dimensions();
  auto dparams = std::vector<Q>{};
  auto dactions = std::vector<Q>{};
  {
    auto qzero = handlers.qq().copy(qparams.front());
    handlers.qq().fill(0, qzero);
    for (size_t i = 0; i < nD; ++i) {
      dparams.emplace_back(handlers.qq().copy(qzero));
      dactions.emplace_back(handlers.qq().copy(qzero));
    }
  }
  const auto oP = 0;
  const auto oQ = dims.nP;
  const auto nQ = dims.nQ - nQdelete;
  const auto nR = rparams.size();
  const auto oR = oQ + nQ;
  const auto oQdelete = oR + nR;
  const auto oDold = oQdelete + nQdelete;
  auto q_indices_new = std::vector<int>(nQ);
  for (size_t j = 0, k = 0; j < dims.nQ; ++j) {
    if (std::find(begin(q_indices_remove), end(q_indices_remove), j) == end(q_indices_remove))
      q_indices_new[k++] = j;
  }
  for (size_t i = 0; i < nD; ++i) {
    for (size_t j = 0; j < dims.nP; ++j) {
      handlers.qp().axpy(lin_trans(i, oP + j), pparams.at(j), dparams[i]);
      handlers.qp().axpy(lin_trans(i, oP + j), pactions.at(j), dactions[i]);
    }
    for (size_t j = 0; j < nQ; ++j) {
      const auto jj = q_indices_new[j];
      handlers.qq().axpy(lin_trans(i, oQ + j), qparams.at(jj), dparams[i]);
      handlers.qq().axpy(lin_trans(i, oQ + j), qactions.at(jj), dactions[i]);
    }
    for (size_t j = 0; j < nR; ++j) {
      handlers.qr().axpy(lin_trans(i, oR + j), rparams.at(j), dparams[i]);
    }
    for (size_t j = 0; j < nQdelete; ++j) {
      const auto jj = q_indices_remove[j];
      handlers.qq().axpy(lin_trans(i, oQdelete + j), qparams.at(jj), dparams[i]);
      handlers.qq().axpy(lin_trans(i, oQdelete + j), qactions.at(jj), dactions[i]);
    }
    for (size_t j = 0; j < dims.nD; ++j) {
      handlers.qq().axpy(lin_trans(i, oDold + j), dparams_old.at(j), dparams[i]);
      handlers.qq().axpy(lin_trans(i, oDold + j), dactions_old.at(j), dactions[i]);
    }
  }
  auto lin_trans_only_R = subspace::Matrix<value_type>({nD, nR});
  lin_trans_only_R.slice() = lin_trans.slice({0, oR}, {nD, oR + nR});
  return std::make_tuple(std::move(dparams), std::move(dactions), lin_trans_only_R);
}

/*!
 * @brief Normalises D parameters and remove any that are null
 */
template <class Q, typename value_type, typename value_type_abs>
void normalise_and_remove_null_D_params(std::vector<Q>& dparams, std::vector<Q>& dactions,
                                        subspace::Matrix<value_type>& lin_trans_D_only_R,
                                        subspace::Matrix<value_type>& solutions_D, const value_type_abs norm_thresh,
                                        array::ArrayHandler<Q, Q>& handler, Logger& logger) {
  logger.msg("normalise_and_remove_null_D_params()", Logger::Trace);
  const auto nD = dparams.size();
  auto to_remove = std::vector<size_t>{};
  for (size_t i = 0; i < nD; ++i) {
    auto norm = handler.dot(dparams[i], dparams[i]);
    norm = std::sqrt(std::abs(norm));
    if (norm > norm_thresh) {
      handler.scal(1. / norm, dparams[i]);
      handler.scal(1. / norm, dactions[i]);
      lin_trans_D_only_R.row(i).scal(1. / norm);
    } else {
      to_remove.push_back(i);
      std::stringstream ss;
      ss << std::setprecision(3) << "remove D parameter i = " << i << ", norm = " << norm;
      logger.msg(ss.str(), Logger::Debug);
    }
  }
  std::sort(to_remove.begin(), to_remove.end(), std::greater());
  for (const auto i : to_remove) {
    dparams.erase(dparams.begin() + i);
    dactions.erase(dactions.begin() + i);
    lin_trans_D_only_R.remove_row(i);
    solutions_D.remove_row(i);
  }
}

/*!
 * @brief Constructs overlap of P+Q+R+Dnew subspace by extending P+Q+R overlap with components from Dnew
 * @param overlap_PQDR overlap of P+Q+D+R subspace where D is the old D space
 * @param dims dimensions of the old subspace which specify P+Q+D distribution
 * @param pparams P space parameters
 * @param qparams Q space parameters
 * @param rparams new R space parameters
 * @param dparams new D space parameters
 * @param handlers array handlers
 * @param logger logger
 * @return overlap of P+Q+R+Dnew subspace
 */
template <class R, class Q, class P, typename value_type>
auto construct_overlap_of_new_subspace(const subspace::Matrix<value_type>& overlap_PQDR,
                                       const subspace::Dimensions& dims, const CVecRef<P> pparams,
                                       const CVecRef<Q>& qparams, const CVecRef<R>& rparams, const CVecRef<Q>& dparams,
                                       ArrayHandlers<R, Q, P>& handlers, Logger& logger) {
  logger.msg("construct_overlap_of_new_subspace()", Logger::Trace);
  const size_t oP = 0, nP = dims.nP;
  const size_t oQ = nP, nQ = dims.nQ;
  const auto oR = nP + nQ, nR = rparams.size();
  const auto oD = oR + nR, nD = dparams.size();
  const auto nX = nP + nQ + nR + nD;
  auto ov = subspace::Matrix<value_type>({nX, nX});
  ov.slice({oP, oP}, {oR, oR}) = overlap_PQDR.slice({oP, oP}, {oR, oR});
  ov.slice({oP, oR}, {oR, oR + nR}) = overlap_PQDR.slice({oP, oR + dims.nD}, {oR, oR + dims.nD + nR});
  ov.slice({oR, oR}, {oR + nR, oR + nR}) =
      overlap_PQDR.slice({oR + dims.nD, oR + dims.nD}, {oR + dims.nD + nR, oR + dims.nD + nR});
  ov.slice({oP, oD}, {oP + nP, oD + nD}) = subspace::util::overlap(pparams, dparams, handlers.qp());
  ov.slice({oQ, oD}, {oQ + nQ, oD + nD}) = subspace::util::overlap(qparams, dparams, handlers.qq());
  ov.slice({oR, oD}, {oR + nR, oD + nD}) = subspace::util::overlap(rparams, dparams, handlers.rq());
  ov.slice({oD, oD}, {oD + nD, oD + nD}) = subspace::util::overlap(dparams, handlers.qq());
  for (size_t i = 0; i < nX; ++i)
    for (size_t j = 0; j < i; ++j)
      ov(i, j) = ov(j, i);
  if (logger.data_dump) {
    logger.msg("overlap P+Q+R+Dnew = " + as_string(ov), Logger::Trace);
  }
  return ov;
}

/*!
 * @brief Apply SVD to check that the subspace is not overcomplete due to D space. If it is mark more Q parameters for
 * removal.
 *
 * For overcomplete problems the D space can be partially null. This routine aims to discover the null space and
 * eliminate it by removing some Q parameters so that a stable D space can be reconstructed.
 *
 * This is done by looking at the smallest right singular vectors and finding D parameters that have the largest
 * weight. The corresponding solution is than searched for the largest contributing Q parameter, which is marked
 * for deletion. If there are no Q parameters, than the D parameter itself is marked for removal
 *
 * @param overlap_PQRD overlap of the subspace, where Q includes parameters marked for removal
 * @param solutions_D solution matrix corresponding to the D parameters
 * @param q_indices_remove Q parameters that are marked for removal so far
 * @param oQ offset to the start of Q block
 * @param nQ total number of Q parameters
 * @return whether the subspace is well conditioned (all singular values are above threshold), and which D parameters to
 * remove
 */
template <typename value_type, typename value_type_abs>
auto fix_overcompleteness(subspace::Matrix<value_type> overlap_PQRD, const subspace::Matrix<value_type>& solutions_D,
                          std::vector<int>& q_indices_remove, size_t oQ, size_t nQ, value_type_abs svd_thresh,
                          Logger& logger) {
  logger.msg("fix_overcompleteness()", Logger::Trace);
  std::sort(begin(q_indices_remove), end(q_indices_remove), std::greater());
  auto q_indices = std::vector<int>{};
  for (size_t i = 0; i < nQ; ++i)
    if (std::find(begin(q_indices_remove), end(q_indices_remove), i) == end(q_indices_remove))
      q_indices.push_back(i);
  for (auto i : q_indices_remove)
    overlap_PQRD.remove_row_col(oQ + i, oQ + i);
  const auto nD = solutions_D.rows();
  const auto oD = overlap_PQRD.cols() - int(nD);
  auto d_indices = std::vector<int>(nD);
  std::iota(begin(d_indices), end(d_indices), 0);
  auto svd_vecs = svd_system(overlap_PQRD.rows(), overlap_PQRD.cols(),
                             array::Span(&overlap_PQRD(0, 0), overlap_PQRD.size()), svd_thresh);
  auto d_indices_remove = std::vector<int>{};
  auto well_conditioned = svd_vecs.empty() || q_indices.empty();
  if (not well_conditioned) {
    for (const auto& svd : svd_vecs) {
      if (!q_indices.empty()) {
        auto contrib = std::vector<value_type_abs>{};
        for (auto i : d_indices)
          contrib.push_back(std::abs(svd.v.at(oD + i)));
        auto itd = std::max_element(begin(contrib), end(contrib));
        auto id = std::distance(begin(contrib), itd);
        auto sol_contrib = std::vector<value_type_abs>{};
        for (auto i : q_indices)
          sol_contrib.push_back(std::abs(solutions_D(d_indices.at(id), oQ + i)));
        auto itq = std::max_element(begin(sol_contrib), end(sol_contrib));
        auto iq = std::distance(begin(sol_contrib), itq);
        {
          auto ss = std::stringstream{};
          ss << "d parameter forming the null space = " << d_indices.at(id) << ", svd value = " << std::setprecision(3)
             << svd.value << ", svd vector contribution = " << *itd;
          logger.msg(ss.str(), Logger::Debug);
          ss = std::stringstream{};
          ss << std::setprecision(3) << "q parameter with max contrib to solution = " << q_indices.at(iq)
             << ", contribution = " << *itq;
          logger.msg(ss.str(), Logger::Debug);
        }
        q_indices_remove.push_back(q_indices.at(iq));
        q_indices.erase(begin(q_indices) + iq);
        d_indices.erase(begin(d_indices) + id);
      }
    }
  } else {
    for (const auto& svd : svd_vecs) {
      auto id = d_indices.back();
      auto ss = std::stringstream{};
      ss << std::setprecision(3) << "removed d parameter, index = " << std::to_string(d_indices.at(id))
         << ", svd.v[i] = " << svd.v.at(oD + d_indices.at(id)) << ", svd.v = ";
      logger.msg(ss.str(), std::begin(svd.v), std::end(svd.v), Logger::Debug);
      d_indices_remove.push_back(d_indices.at(id));
      d_indices.erase(begin(d_indices) + id);
    }
  }
  return std::make_tuple(well_conditioned, d_indices_remove);
}

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
auto get_new_working_set(const std::vector<int>& working_set, const CVecRef<R>& params, const CVecRef<R>& wparams) {
  auto new_indices = find_ref(wparams, params);
  auto new_working_set = std::vector<int>{};
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
auto propose_rspace(LinearEigensystem<R, Q, P>& solver, const VecRef<R>& parameters, const VecRef<R>& residuals,
                    subspace::XSpaceI<R, Q, P>& xspace, const subspace::Matrix<value_type>& solutions,
                    ArrayHandlers<R, Q, P>& handlers, Logger& logger, value_type_abs svd_thresh,
                    value_type_abs res_norm_thresh, int max_size_qspace) {
  logger.msg("itsolv::detail::propose_rspace", Logger::Trace);
  logger.msg("dimensions {nP, nQ, nD, nW} = " + std::to_string(xspace.dimensions().nP) + ", " +
                 std::to_string(xspace.dimensions().nQ) + ", " + std::to_string(xspace.dimensions().nD) + ", " +
                 std::to_string(solver.working_set().size()),
             Logger::Trace);
  auto wresidual = wrap<R>(residuals.begin(), residuals.begin() + solver.working_set().size());
  normalise(wresidual, handlers.rr(), logger);
  auto null_param_indices = subspace::util::modified_gram_schmidt(wresidual, handlers.rr(), 1.0e-14);
  for (auto it = null_param_indices.rbegin(); it != null_param_indices.rend(); ++it)
    wresidual.erase(begin(wresidual) + (*it));
  // propose working space by orthogonalising against P+Q
  auto overlap_PQDR = append_overlap_with_r(xspace.data.at(subspace::EqnData::S), cwrap(wresidual), xspace.cparamsp(),
                                            xspace.cparamsq(), handlers, logger);
  auto [q_indices_remove, lin_trans, norm] = calculate_transformation_to_orthogonal_rspace(
      overlap_PQDR, solutions, xspace.dimensions(), logger, svd_thresh, res_norm_thresh, max_size_qspace);
  if (logger.data_dump) {
    logger.msg("overlap P+Q+Z = " + subspace::as_string(overlap_PQDR), Logger::Info);
    logger.msg("linear transformation = " + subspace::as_string(lin_trans), Logger::Info);
    logger.msg("norm = ", norm.begin(), norm.end(), Logger::Debug);
  }
  logger.msg("remove Q space indices = ", q_indices_remove.begin(), q_indices_remove.end(), Logger::Debug);
  auto wparams = wrap<R>(parameters.begin(), parameters.begin() + wresidual.size());
  auto qparams_new = remove_elements(xspace.cparamsq(), q_indices_remove);
  construct_orthonormal_Rparams(wparams, wresidual, lin_trans, norm, xspace.cparamsp(), qparams_new, handlers);
  normalise(wparams, handlers.rr(), logger);
  auto params_qd = xspace.cparamsq();
  auto paramsd = xspace.cparamsd();
  std::copy(begin(paramsd), end(paramsd), std::back_inserter(params_qd));
  overlap_PQDR = append_overlap_with_r(xspace.data.at(subspace::EqnData::S), cwrap(wparams), xspace.cparamsp(),
                                       params_qd, handlers, logger);
  std::vector<Q> dparams, dactions;
  subspace::Matrix<value_type> lin_trans_D_only_R;
  for (bool well_conditioned = false; not well_conditioned;) {
    logger.msg("proposing D space", Logger::Debug);
    auto [lin_trans_D, solutions_D] = dspace::propose_dspace(solutions, xspace.dimensions(), q_indices_remove,
                                                             overlap_PQDR, wparams.size(), res_norm_thresh, logger);
    std::tie(dparams, dactions, lin_trans_D_only_R) =
        construct_orthonormal_Dparams(xspace, lin_trans_D, q_indices_remove, cwrap(wparams), handlers, logger);
    normalise_and_remove_null_D_params(dparams, dactions, lin_trans_D_only_R, solutions_D, res_norm_thresh,
                                       handlers.qq(), logger);
    auto overlap_PQRDnew =
        construct_overlap_of_new_subspace(overlap_PQDR, xspace.dimensions(), xspace.cparamsp(), xspace.cparamsq(),
                                          cwrap(wparams), cwrap(dparams), handlers, logger);
    if (logger.data_dump) {
      logger.msg("D params in subspace = " + subspace::as_string(lin_trans_D), Logger::Info);
    }
    std::vector<int> d_indices_remove;
    std::tie(well_conditioned, d_indices_remove) =
        fix_overcompleteness(overlap_PQRDnew, solutions_D, q_indices_remove, xspace.dimensions().oQ,
                             xspace.dimensions().nQ, svd_thresh, logger);
    std::sort(begin(d_indices_remove), end(d_indices_remove), std::greater());
    for (auto id : d_indices_remove) {
      dparams.erase(dparams.begin() + id);
      dactions.erase(dactions.begin() + id);
      lin_trans_D_only_R.remove_row(id);
    }
    if (!d_indices_remove.empty())
      logger.msg("remove d indices ", std::begin(d_indices_remove), std::end(d_indices_remove), Logger::Info);
    logger.msg("subspace is overcomplete = " + std::to_string(!well_conditioned), Logger::Debug);
  }
  std::sort(begin(q_indices_remove), end(q_indices_remove), std::greater());
  for (auto iq : q_indices_remove)
    xspace.eraseq(iq);
  auto wdparams = wrap(dparams);
  auto wdactions = wrap(dactions);
  xspace.update_dspace(wdparams, wdactions, lin_trans_D_only_R);
  auto new_working_set = get_new_working_set(solver.working_set(), cwrap(residuals), cwrap(wresidual));
  return new_working_set;
}
} // namespace molpro::linalg::itsolv::detail

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_PROPOSE_RSPACE_H
