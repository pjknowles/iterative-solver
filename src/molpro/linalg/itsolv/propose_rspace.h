#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_PROPOSE_RSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_PROPOSE_RSPACE_H
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/helper.h>
#include <molpro/linalg/itsolv/subspace/Dimensions.h>
#include <molpro/linalg/itsolv/subspace/ISubspaceSolver.h>
#include <molpro/linalg/itsolv/subspace/IXSpace.h>
#include <molpro/linalg/itsolv/subspace/QSpace.h>
#include <molpro/linalg/itsolv/subspace/gram_schmidt.h>
#include <molpro/linalg/itsolv/subspace/util.h>
#include <molpro/linalg/itsolv/util.h>
#include <molpro/linalg/itsolv/wrap_util.h>
#include <molpro/profiler/Profiler.h>

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
          overlap_proj(i, ii) +=
              solutions_proj(i, nQd + j) * solutions_proj(ii, k) * overlap(dims.oD + j, dims.oQ + remove_qspace[k]);
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
      std::stringstream ss;
      ss << std::setprecision(3) << "remove projected solution parameter i = " << i << ", norm = " << norm_proj[i];
      logger.msg(ss.str(), Logger::Info);
    }
  }
  if (logger.data_dump) {
    logger.msg("norm  = ", std::begin(norm_proj), std::end(norm_proj), Logger::Info);
    logger.msg("parameters after normalisation = " + as_string(parameters), Logger::Info);
    logger.msg("overlap of parameters = " + as_string(overlap), Logger::Info);
  }
}

/*!
 * @brief Transforms to a stable subspace of projected solutions via SVD
 *
 * Performs SVD and transforms to right-singular vectors with singular values greater than threshold.
 *
 * @param solutions_proj solutions projected onto Qd+D space
 * @param solutions solutions in the current subspace
 * @param overlap_proj overlap matrix of projected solutions
 * @param norm_thresh
 * @param logger
 * @returns stable subspace of projected solutions
 */
template <typename value_type, typename value_type_abs>
auto remove_null_projected_solutions(const subspace::Matrix<value_type>& solutions_proj,
                                     const subspace::Matrix<value_type>& overlap_proj, const value_type_abs svd_thresh,
                                     Logger& logger) {
  logger.msg("remove_null_projected_solutions()", Logger::Trace);
  logger.msg("nS on entry = " + std::to_string(solutions_proj.rows()), Logger::Debug);
  value_type* m = const_cast<std::vector<value_type>&>(overlap_proj.data()).data();
  auto svd_vecs = svd_system(overlap_proj.rows(), overlap_proj.cols(), array::Span(m, overlap_proj.size()),
                             std::numeric_limits<value_type_abs>::max());
  svd_vecs.remove_if([&svd_thresh](const auto& el) { return el.value < svd_thresh; });
  svd_vecs.sort([](const auto& lt, const auto& rt) { return lt.value < rt.value; });
  const auto nD = svd_vecs.size();
  const auto nX = solutions_proj.cols();
  auto solutions_stable = subspace::Matrix<value_type>({nD, nX});
  auto svd = svd_vecs.begin();
  for (size_t i = 0; i < nD; ++i, ++svd)
    for (size_t j = 0; j < overlap_proj.cols(); ++j)
      for (size_t k = 0; k < nX; ++k)
        solutions_stable(i, k) += svd->v[j] * solutions_proj(j, k);
  logger.msg("nS without null space = " + std::to_string(solutions_stable.rows()), Logger::Debug);
  return solutions_stable;
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
auto construct_full_subspace_overlap(const subspace::Matrix<value_type>& solutions_proj,
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
} // namespace dspace

/*!
 * @brief Constructs overlap of the full subspace by appending overlap with new parameters to the overlap of previous
 * subspace
 *
 * @param ov overlap of P+Q+D
 * @param params new parameters
 * @param pparams P space parameters
 * @param qparams Q space parameters
 * @param dparams D space parameters
 * @param handlers array handlers
 * @param logger logger
 */
template <class R, class Q, class P, typename value_type>
auto append_overlap_with_r(const subspace::Matrix<value_type>& overlap, const CVecRef<R>& params,
                           const CVecRef<P>& pparams, const CVecRef<Q>& qparams, const CVecRef<Q>& dparams,
                           ArrayHandlers<R, Q, P>& handlers, Logger& logger) {
  logger.msg("append_overlap_with_r()", Logger::Trace);
  const auto nP = pparams.size(), nQ = qparams.size(), nD = dparams.size(), nN = params.size();
  const auto nX = nP + nQ + nD + nN;
  const auto oP = 0;
  const auto oQ = oP + nP;
  const auto oD = oQ + nQ;
  const auto oN = oD + nD;
  auto ov = overlap;
  ov.resize({nX, nX}); // solutions are last
  ov.slice({oN, oN}, {oN + nN, oN + nN}) = subspace::util::overlap(params, handlers.rr());
  ov.slice({oN, oP}, {oN + nN, oP + nP}) = subspace::util::overlap(params, pparams, handlers.rp());
  ov.slice({oN, oQ}, {oN + nN, oQ + nQ}) = subspace::util::overlap(params, qparams, handlers.rq());
  ov.slice({oN, oD}, {oN + nN, oD + nD}) = subspace::util::overlap(params, dparams, handlers.rq());
  auto copy_upper_to_lower = [&ov, oN, nN](size_t oX, size_t nX) {
    for (size_t i = 0; i < nX; ++i)
      for (size_t j = 0; j < nN; ++j)
        ov(oX + i, oN + j) = ov(oN + j, oX + i);
  };
  copy_upper_to_lower(oP, nP);
  copy_upper_to_lower(oQ, nQ);
  copy_upper_to_lower(oD, nD);
  if (logger.data_dump) {
    logger.msg("full overlap P+Q+D+R = " + as_string(ov), Logger::Info);
  }
  return ov;
}

/*!
 * @brief Ensures that size Q space is within limit by proposing Q parameters for deletion
 * @param dims dimensions of current P+Q+D subspace
 * @param max_size_qspace maximum size of Q space
 * @param solutions solutions matrix
 * @param logger logger
 * @return q parameters marked for deletions
 */
template <typename value_type>
auto limit_qspace_size(const subspace::Dimensions& dims, const size_t max_size_qspace,
                       const subspace::Matrix<value_type>& solutions, Logger& logger) {
  logger.msg("limit_qspace_size()", Logger::Trace);
  using value_type_abs = decltype(std::abs(solutions(0, 0)));
  auto q_delete = std::vector<int>{};
  auto q_indices = std::vector<int>(dims.nQ);
  std::iota(begin(q_indices), end(q_indices), 0);
  const auto nSol = solutions.rows();
  while (q_indices.size() > max_size_qspace) {
    auto max_contrib_to_solution = std::vector<value_type_abs>{};
    for (auto i : q_indices) {
      auto contrib = std::vector<value_type_abs>(nSol);
      for (size_t j = 0; j < nSol; ++j)
        contrib[j] = std::abs(solutions(j, dims.oQ + i));
      max_contrib_to_solution.emplace_back(*std::max_element(begin(contrib), end(contrib)));
    }
    auto it_min = std::min_element(begin(max_contrib_to_solution), end(max_contrib_to_solution));
    auto i = std::distance(begin(max_contrib_to_solution), it_min);
    q_delete.push_back(q_indices.at(i));
    q_indices.erase(begin(q_indices) + i);
    logger.msg("contribution to solutions =", std::begin(max_contrib_to_solution), std::end(max_contrib_to_solution),
               Logger::Info);
    logger.msg("delete Q i = " + std::to_string(i), Logger::Info);
  }
  return q_delete;
}

/*!
 * @brief Constructs the new D space by projecting the solutions onto Qd+D subspace and ensuring they are well
 * conditioned.
 * @param xspace X space container
 * @param q_delete Q parameters marked for deletion
 * @param norm_thresh norm threshold for parameters to be considered null
 * @param svd_thresh svd threshold for marking the null space
 * @param handler array handler
 * @param logger logger
 * @return D space parameters and actions
 */
template <class R, class Q, class P, typename value_type, typename value_type_abs>
auto construct_dspace(const subspace::Matrix<value_type>& solutions, const subspace::IXSpace<R, Q, P>& xspace,
                      const std::vector<int>& q_delete, const value_type_abs norm_thresh,
                      const value_type_abs svd_thresh, array::ArrayHandler<Q, Q>& handler, Logger& logger) {
  const auto dims = xspace.dimensions();
  const auto overlap = xspace.data.at(subspace::EqnData::S);
  auto solutions_proj = dspace::construct_projected_solution(solutions, dims, q_delete, logger);
  auto overlap_proj = dspace::construct_projected_solutions_overlap(solutions_proj, overlap, dims, q_delete, logger);
  dspace::remove_null_norm_and_normalise(solutions_proj, overlap_proj, norm_thresh, logger);
  solutions_proj = dspace::remove_null_projected_solutions(solutions_proj, overlap_proj, svd_thresh, logger);
  overlap_proj = dspace::construct_projected_solutions_overlap(solutions_proj, overlap, dims, q_delete, logger);
  dspace::remove_null_norm_and_normalise(solutions_proj, overlap_proj, norm_thresh, logger);
  auto overlap_full_subspace = dspace::construct_full_subspace_overlap(solutions_proj, dims, q_delete, overlap, 0);
  auto svd_vecs = svd_system(overlap_full_subspace.rows(), overlap_full_subspace.cols(),
                             array::Span(&overlap_full_subspace(0, 0), overlap_full_subspace.size()), svd_thresh);
  assert(svd_vecs.empty() && "P+Q+D subspace should be stable by construction");
  const auto nD = solutions_proj.rows();
  const auto nQd = q_delete.size();
  assert(nQd + dims.nD == solutions_proj.cols());
  const auto &qparams = xspace.cparamsq(), qactions = xspace.cactionsq();
  const auto &dparams = xspace.cparamsd(), dactions = xspace.cactionsd();
  std::vector<Q> dparams_new, dactions_new;
  {
    // FIXME need a simpler mechanism for constructing null parameters
    Q const* q = nullptr;
    if (!qparams.empty())
      q = &qparams.front().get();
    else if (!dparams.empty())
      q = &dparams.front().get();
    if (q) {
      for (size_t i = 0; i < nD; ++i) {
        dparams_new.emplace_back(handler.copy(*q));
        dactions_new.emplace_back(handler.copy(*q));
        handler.fill(0, dparams_new.back());
        handler.fill(0, dactions_new.back());
      }
    }
  }
  for (size_t i = 0; i < nD; ++i) {
    for (size_t j = 0; j < q_delete.size(); ++j) {
      handler.axpy(solutions_proj(i, j), qparams.at(q_delete[j]), dparams_new.at(i));
      handler.axpy(solutions_proj(i, j), qactions.at(q_delete[j]), dactions_new.at(i));
    }
    for (size_t j = 0; j < dims.nD; ++j) {
      handler.axpy(solutions_proj(i, nQd + j), dparams.at(j), dparams_new.at(i));
      handler.axpy(solutions_proj(i, nQd + j), dactions.at(j), dactions_new.at(i));
    }
  }
  for (size_t i = 0; i < nD; ++i) {
    auto norm = std::sqrt(std::abs(handler.dot(dparams_new.at(i), dparams_new.at(i))));
    handler.scal(1. / norm, dparams_new[i]);
    handler.scal(1. / norm, dactions_new[i]);
  }
  return std::make_tuple(std::move(dparams_new), std::move(dactions_new));
}

/*!
 * @brief Orthogonalises R parameters against P+Q+D subspace (and themselves)
 *
 * New vectors with norm less than threshold are considered null and are not normalised.
 * Their indices in params are returned.
 *
 * @param rparams R space parameters
 * @param overlap overlap of P+Q+D subspace
 * @param dims dimensions of P+Q+D subspace
 * @param pparams P space parameters
 * @param qparams Q space parameters
 * @param dparams D space parameters
 * @param handlers array handlers
 * @param norm_thresh  parameters with norm less than threshold are considered null and won't be orthogonalised against
 * @return indices of null parameters
 */
template <class R, class Q, class P, typename value_type, typename value_type_abs>
auto modified_gram_schmidt(const VecRef<R>& rparams, const subspace::Matrix<value_type>& overlap,
                           const subspace::Dimensions& dims, const CVecRef<P>& pparams, const CVecRef<Q>& qparams,
                           const CVecRef<Q>& dparams, const value_type_abs norm_thresh,
                           ArrayHandlers<R, Q, P>& handlers, Logger& logger) {
  logger.msg("modified_gram_schmidt()", Logger::Trace);
  const auto nR = rparams.size(), nP = pparams.size(), nQ = qparams.size(), nD = dparams.size();
  // std::cout << "modified_gram_schmidt() nR="<<nR<<" nQ="<<nQ<<" nD="<<nD<<std::endl;
  assert(nP == dims.nP && nQ == dims.nQ && nD == dims.nD);
  auto orthogonalise = [&overlap, &rparams, nR](const auto& xparams, auto& handler, const size_t oX, const size_t nX) {
    for (size_t i = 0; i < nX; ++i) {
      auto norm = std::abs(overlap(oX + i, oX + i));
      if (nR > 0) {
        auto dot_mat = handler.gemm_inner(cwrap(rparams), cwrap_arg(xparams.at(i).get()));
        std::pair<size_t, size_t> mcoeff_dim = std::make_pair(1, nR);
        Matrix<double> mcoeff(dot_mat.data(), mcoeff_dim);
        for (size_t j = 0; j < nR; ++j) {
          mcoeff(0, j) = -mcoeff(0, j) / norm;
        }
        handler.gemm_outer(mcoeff, cwrap_arg(xparams.at(i).get()), rparams);
      }
    }
  };
  auto prof = molpro::Profiler::single();
  prof->start("orthoganalise");
  orthogonalise(pparams, handlers.rp(), dims.oP, nP);
  orthogonalise(qparams, handlers.rq(), dims.oQ, nQ);
  orthogonalise(dparams, handlers.rq(), dims.oD, nD);
  prof->stop();
  prof->start("get null_params");
  auto null_params = std::vector<int>{};
  for (size_t i = 0; i < nR; ++i) {
    auto norm = std::sqrt(std::abs(handlers.rr().dot(rparams[i], rparams[i])));
    if (norm > norm_thresh) {
      handlers.rr().scal(1. / norm, rparams[i]);
      for (size_t j = i + 1; j < nR; ++j) {
        auto ov = handlers.rr().dot(rparams[i], rparams[j]);
        handlers.rr().axpy(-ov, rparams[i], rparams[j]);
      }
    } else {
      null_params.push_back(i);
    }
  }
  prof->stop();
  return null_params;
}

/*!
 * @brief Deduces a set of parameters that are redundant due to linear dependencies
 *
 * Only the last nR parameters are considered for removal. Linea dependencies are discovered by performing SVD of the
 * overlap matrix.
 *
 * @param overlap overlap matrix of the full subspace
 * @param oR offset to the start of parameter block
 * @param nR number of parameters to consider for removal
 * @param svd_thresh singular value threshold for choosing the null space
 * @param logger logger
 * @return indices of the last nR parameters that are considered redundant.
 */
template <typename value_type, typename value_type_abs>
auto redundant_parameters(const subspace::Matrix<value_type>& overlap, const size_t oR, const size_t nR,
                          const value_type_abs svd_thresh, Logger& logger) {
  auto prof = molpro::Profiler::single();
  prof->start("itsolv::svd_system");
  logger.msg("redundant_parameters()", Logger::Trace);
  auto redundant_params = std::vector<int>{};
  auto rspace_indices = std::vector<int>(nR);
  std::iota(std::begin(rspace_indices), std::end(rspace_indices), 0);
  auto svd = svd_system(overlap.rows(), overlap.cols(),
                        array::Span(const_cast<value_type*>(overlap.data().data()), overlap.size()), svd_thresh);
  prof->stop();
  prof->start("find redundant parameters");
  for (const auto& singular_system : svd) {
    if (!rspace_indices.empty()) {
      auto rspace_contribution = std::vector<value_type_abs>{};
      for (auto i : rspace_indices)
        rspace_contribution.push_back(std::abs(singular_system.v.at(oR + i)));
      auto it_min = std::max_element(std::begin(rspace_contribution), std::end(rspace_contribution));
      auto imin = std::distance(std::begin(rspace_contribution), it_min);
      redundant_params.push_back(rspace_indices[imin]);
      rspace_indices.erase(std::begin(rspace_indices) + imin);
      std::stringstream ss;
      ss << std::setprecision(3) << "redundant parameter found, i = " << redundant_params.back()
         << ", svd.value = " << singular_system.value
         << ", svd.v[i] = " << singular_system.v[oR + redundant_params.back()];
      logger.msg(ss.str(), Logger::Info);
    }
  }
  prof->stop();
  return redundant_params;
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
template <class R, class Q, class P, typename value_type_abs>
auto propose_rspace(IterativeSolver<R, Q, P>& solver, const VecRef<R>& parameters, const VecRef<R>& residuals,
                    subspace::IXSpace<R, Q, P>& xspace, subspace::ISubspaceSolver<R, Q, P>& subspace_solver,
                    ArrayHandlers<R, Q, P>& handlers, Logger& logger, value_type_abs svd_thresh,
                    value_type_abs res_norm_thresh, int max_size_qspace, molpro::profiler::Profiler& profiler) {
  // auto prof = profiler.push("itsolv::propose_rspace"); // FIXME two separate profilers
  auto prof = molpro::Profiler::single();
  logger.msg("itsolv::detail::propose_rspace", Logger::Trace);
  logger.msg("dimensions {nP, nQ, nD, nW} = " + std::to_string(xspace.dimensions().nP) + ", " +
                 std::to_string(xspace.dimensions().nQ) + ", " + std::to_string(xspace.dimensions().nD) + ", " +
                 std::to_string(solver.working_set().size()),
             Logger::Trace);
  profiler.start("itsolv::ISubspaceSolver::solutions");
  auto solutions = subspace_solver.solutions();
  profiler.stop();
  auto q_delete = limit_qspace_size(xspace.dimensions(), max_size_qspace, solutions, logger);
  logger.msg("delete Q parameter indices = ", q_delete.begin(), q_delete.end(), Logger::Debug);
  if (!q_delete.empty()) {
    auto prof = profiler.push("construct_dspace");
    auto [dparams, dactions] =
        construct_dspace(solutions, xspace, q_delete, res_norm_thresh, svd_thresh, handlers.qq(), logger);
    std::sort(begin(q_delete), end(q_delete), std::greater<int>());
    for (auto iq : q_delete)
      xspace.eraseq(iq);
    auto wdparams = wrap(dparams);
    auto wdactions = wrap(dactions);
    xspace.update_dspace(wdparams, wdactions);
    auto eigenvalues_ref = subspace_solver.eigenvalues();
    subspace_solver.solve(xspace, solutions.rows());
    auto eigval_error = std::vector<double>{};
    std::transform(std::begin(eigenvalues_ref), std::end(eigenvalues_ref), std::begin(subspace_solver.eigenvalues()),
                   std::back_inserter(eigval_error), [](auto& e_ref, auto& e_new) { return std::abs(e_ref - e_new); });
    logger.msg("eigenvalue error due to new D space = ", std::begin(eigval_error), std::end(eigval_error),
               Logger::Debug);
    // FIXME Optionally, solve the subspace problem again and get an estimate of the error due to new D
  }
  // Use modified GS to orthonormalise z against P+Q+D, removing any null parameters.
  auto wresidual = wrap(residuals.begin(), residuals.begin() + solver.working_set().size());
  profiler.start("normalise");
  normalise(wresidual, handlers.rr(), logger);
  profiler.stop();
  profiler.start("append_overlap_with_r");
  prof->start("append_overlap_with_r");
  const auto full_overlap =
      append_overlap_with_r(xspace.data.at(subspace::EqnData::S), cwrap(wresidual), xspace.cparamsp(),
                            xspace.cparamsq(), xspace.cparamsd(), handlers, logger);
  profiler.stop();
  prof->stop();
  prof->start("redundant_indices");
  auto redundant_indices =
      redundant_parameters(full_overlap, xspace.dimensions().nX, wresidual.size(), svd_thresh, logger);
  prof->stop();
  logger.msg("redundant indices = ", std::begin(redundant_indices), std::end(redundant_indices), Logger::Debug);
  util::delete_parameters(redundant_indices, wresidual);
  profiler.start("modified_gram_schmidt");
  prof->start("modified_gram_schmidt");
  auto null_param_indices =
      modified_gram_schmidt(wresidual, xspace.data.at(subspace::EqnData::S), xspace.dimensions(), xspace.cparamsp(),
                            xspace.cparamsq(), xspace.cparamsd(), res_norm_thresh, handlers, logger);
  profiler.stop();
  prof->stop();
  // Now that there is SVD null_param_indices should always be empty
  logger.msg("null parameters = ", std::begin(null_param_indices), std::end(null_param_indices), Logger::Debug);
  util::delete_parameters(null_param_indices, wresidual);
  normalise(wresidual, handlers.rr(), logger);
  for (size_t i = 0; i < wresidual.size(); ++i)
    handlers.rr().copy(parameters.at(i), wresidual.at(i));
  profiler.start("get_new_working_set");
  auto new_working_set = get_new_working_set(solver.working_set(), cwrap(residuals), cwrap(wresidual));
  profiler.stop();
  return new_working_set;
}
} // namespace molpro::linalg::itsolv::detail

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_PROPOSE_RSPACE_H
