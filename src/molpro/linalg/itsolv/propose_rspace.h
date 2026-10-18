#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_PROPOSE_RSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_PROPOSE_RSPACE_H
#include <Eigen/Dense>
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/helper.h>
#include <molpro/linalg/itsolv/qspace_options.h>
#include <molpro/linalg/itsolv/rspace_options.h>
#include <molpro/linalg/itsolv/subspace/Dimensions.h>
#include <molpro/linalg/itsolv/subspace/ISubspaceSolver.h>
#include <molpro/linalg/itsolv/subspace/IXSpace.h>
#include <molpro/linalg/itsolv/subspace/QSpace.h>
#include <molpro/linalg/itsolv/subspace/gram_schmidt.h>
#include <molpro/linalg/itsolv/subspace/util.h>
#include <molpro/linalg/itsolv/util.h>
#include <molpro/linalg/itsolv/wrap_util.h>
#include <molpro/profiler/Profiler.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <format>
#include <functional>
#include <limits>
#include <ranges>
#include <utility>
#include <vector>

namespace molpro::linalg::itsolv::detail {

// Normalises each parameter to unit length, leaving alone only vectors whose length is zero or too small to invert.
// The threshold must not depend on the scale of the problem: a preconditioned residual can legitimately be very short
// (for example with a large diagonal), and whether it is redundant is decided afterwards from the overlap of
// normalised vectors.
template <class R>
void normalise(VecRef<R>& params, array::ArrayHandler<R, R>& handler, Logger& logger,
               typename array::ArrayHandler<R, R>::value_type_abs thresh =
                   1 / std::numeric_limits<typename array::ArrayHandler<R, R>::value_type_abs>::max()) {
  for (auto& p : params) {
    // the length of a vector is a magnitude, real even when the elements are complex
    const typename array::ArrayHandler<R, R>::value_type_abs dot = std::sqrt(std::abs(handler.dot(p, p)));
    if (dot > thresh and std::isfinite(dot)) {
      handler.scal(1. / dot, p);
    } else {
      logger.warn("parameter's length is too small for normalisation, dot = " + std::format("{:.2e}", double(dot)));
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
                                  const std::vector<std::size_t>& remove_qspace, Logger& logger) {
  logger.trace("construct_projected_solution()");
  const auto nQd = remove_qspace.size();
  const auto nSol = solutions.rows();
  auto solutions_proj = subspace::Matrix<value_type>({nSol, nQd + dims.nD});
  for (size_t i = 0; i < nSol; ++i) {
    for (size_t j = 0; j < nQd; ++j) {
      solutions_proj(i, j) = solutions(i, dims.oQ() + remove_qspace[j]);
    }
    for (size_t j = 0; j < dims.nD; ++j) {
      solutions_proj(i, nQd + j) = solutions(i, dims.oD() + j);
    }
  }
  logger.debug("nSol, nQd, nD", nSol, nQd, dims.nD);
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
                                           const subspace::Dimensions& dims,
                                           const std::vector<std::size_t>& remove_qspace, Logger& logger) {
  logger.trace("construct_projected_solution_overlap()");
  const auto nSol = solutions_proj.rows();
  const auto nQd = remove_qspace.size();
  auto overlap_proj = subspace::Matrix<value_type>({nSol, nSol});
  for (size_t i = 0; i < nSol; ++i) {
    for (size_t ii = 0; ii <= i; ++ii) {
      for (size_t j = 0; j < nQd; ++j) {
        for (size_t k = 0; k < nQd; ++k) {
          overlap_proj(i, ii) += solutions_proj(i, j) * solutions_proj(ii, k) *
                                 overlap(dims.oQ() + remove_qspace[j], dims.oQ() + remove_qspace[k]);
        }
        for (size_t k = 0; k < dims.nD; ++k) {
          overlap_proj(i, ii) +=
              solutions_proj(i, j) * solutions_proj(ii, nQd + k) * overlap(dims.oQ() + remove_qspace[j], dims.oD() + k);
        }
      }
      for (size_t j = 0; j < dims.nD; ++j) {
        for (size_t k = 0; k < dims.nD; ++k) {
          overlap_proj(i, ii) +=
              solutions_proj(i, nQd + j) * solutions_proj(ii, nQd + k) * overlap(dims.oD() + j, dims.oD() + k);
        }
        for (size_t k = 0; k < nQd; ++k) {
          overlap_proj(i, ii) +=
              solutions_proj(i, nQd + j) * solutions_proj(ii, k) * overlap(dims.oD() + j, dims.oQ() + remove_qspace[k]);
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
  logger.trace("remove_null_norm_and_normalise()");
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
      logger.info(ss.str());
    }
  }
  logger.data_dump("norm  = ", norm_proj);
  logger.data_dump("parameters after normalisation = ", parameters);
  logger.data_dump("overlap of parameters = ", overlap);
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
  logger.trace("remove_null_projected_solutions()");
  logger.debug("nS on entry = ", solutions_proj.rows());
  value_type* m = const_cast<std::vector<value_type>&>(overlap_proj.data()).data();
  auto svd_vecs = svd_system(overlap_proj.rows(), overlap_proj.cols(), array::Span(m, overlap_proj.size()),
                             std::numeric_limits<value_type_abs>::max(), true);
  svd_vecs.remove_if([&svd_thresh](const auto& el) { return real_part(el.value) < svd_thresh; });
  svd_vecs.sort([](const auto& lt, const auto& rt) { return real_part(lt.value) < real_part(rt.value); });
  const auto nD = svd_vecs.size();
  const auto nX = solutions_proj.cols();
  auto solutions_stable = subspace::Matrix<value_type>({nD, nX});
  auto svd = svd_vecs.begin();
  for (size_t i = 0; i < nD; ++i, ++svd)
    for (size_t j = 0; j < overlap_proj.cols(); ++j)
      for (size_t k = 0; k < nX; ++k)
        solutions_stable(i, k) += svd->v[j] * solutions_proj(j, k);
  logger.debug("nS without null space = ", solutions_stable.rows());
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
                                     const subspace::Dimensions& dims, const std::vector<std::size_t>& remove_qspace,
                                     const subspace::Matrix<value_type>& overlap, const size_t nR) {
  const auto nDnew = solutions_proj.rows();
  const auto nQd = remove_qspace.size();
  const auto nQ = dims.nQ - nQd;
  auto ov = overlap;
  for (size_t i = 0; i < dims.nD; ++i) {
    ov.remove_row_col(dims.oD(), dims.oD());
  }
  auto is_Qdelete = [&remove_qspace](size_t i) {
    return std::find(begin(remove_qspace), end(remove_qspace), i) != end(remove_qspace);
  };
  for (size_t i = 0, j = 0; i < dims.nQ; ++i) {
    if (is_Qdelete(i))
      ov.remove_row_col(dims.oQ() + j, dims.oQ() + j);
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
      ov(oDnew + i, j) += solutions_proj(i, k) * overlap(jj, dims.oQ() + remove_qspace[k]);
    for (size_t k = 0; k < dims.nD; ++k)
      ov(oDnew + i, j) += solutions_proj(i, nQd + k) * overlap(jj, dims.oD() + k);
    ov(j, oDnew + i) = ov(oDnew + i, j);
  };
  for (size_t i = 0; i < nDnew; ++i) {
    for (size_t j = 0; j < dims.nP; ++j)
      accumulate_ov_offdiag(i, j, dims.oP() + j);
    for (size_t j = 0, jj = 0; j < dims.nQ; ++j)
      if (!is_Qdelete(j))
        accumulate_ov_offdiag(i, dims.nP + jj++, dims.oQ() + j);
    for (size_t j = 0; j < nR; ++j)
      accumulate_ov_offdiag(i, dims.nP + nQ + j, dims.nX() + j);
  }
  for (size_t i = 0; i < nDnew; ++i) {
    for (size_t j = 0; j <= i; ++j) {
      for (size_t k = 0; k < nQd; ++k) {
        for (size_t l = 0; l < nQd; ++l)
          ov(oDnew + i, oDnew + j) += solutions_proj(i, k) * solutions_proj(j, l) *
                                      overlap(dims.oQ() + remove_qspace[k], dims.oQ() + remove_qspace[l]);
        for (size_t l = 0; l < dims.nD; ++l)
          ov(oDnew + i, oDnew + j) +=
              solutions_proj(i, k) * solutions_proj(j, nQd + l) * overlap(dims.oQ() + remove_qspace[k], dims.oD() + l);
      }
      for (size_t k = 0; k < dims.nD; ++k) {
        for (size_t l = 0; l < nQd; ++l)
          ov(oDnew + i, oDnew + j) +=
              solutions_proj(i, nQd + k) * solutions_proj(j, l) * overlap(dims.oD() + k, dims.oQ() + remove_qspace[l]);
        for (size_t l = 0; l < dims.nD; ++l)
          ov(oDnew + i, oDnew + j) +=
              solutions_proj(i, nQd + k) * solutions_proj(j, nQd + l) * overlap(dims.oD() + k, dims.oD() + l);
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
  logger.trace("append_overlap_with_r()");
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
  logger.data_dump("full overlap P+Q+D+R = ", ov);
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
std::vector<std::size_t> limit_qspace_size(const subspace::Dimensions& dims, const QSpaceOptions& opts,
                                           const subspace::Matrix<value_type>& solutions, Logger& logger) {
  using value_type_abs = decltype(std::abs(solutions(0, 0)));
  using contrib_pair = std::pair<std::size_t, value_type_abs>;
  using std::begin;
  using std::end;

  logger.trace("limit_qspace_size()");

  assert(opts.max_size >= opts.min_size);

  const bool must_trim_size = dims.nQ > opts.max_size;
  const bool may_trim_size = dims.nQ > opts.min_size && opts.contrib_thresh > 0;

  if (!must_trim_size && !may_trim_size) {
    return {};
  }

  const std::size_t min_deletions = must_trim_size ? dims.nQ - opts.max_size : 0;
  const std::size_t max_deletions = dims.nQ - opts.min_size;

  std::vector<std::size_t> q_delete;
  q_delete.reserve(min_deletions);

  const std::size_t nSol = solutions.rows();

  std::vector<contrib_pair> max_contrib_to_solution;
  max_contrib_to_solution.reserve(nSol);

  for (std::size_t idx : std::ranges::views::iota(std::size_t(0), dims.nQ)) {
    max_contrib_to_solution.emplace_back(std::make_pair(idx, value_type_abs(0)));

    for (size_t j = 0; j < nSol; ++j) {
      const value_type_abs current = std::abs(solutions(j, dims.oQ() + idx));
      max_contrib_to_solution.back().second = std::max(max_contrib_to_solution.back().second, current);
    }
  }

  std::ranges::sort(max_contrib_to_solution, std::less<>{}, &contrib_pair::second);

  logger.trace("contribution to solutions =",
               max_contrib_to_solution | std::ranges::views::transform([](const contrib_pair& p) { return p.second; }));

  for (auto [idx, contrib] : max_contrib_to_solution) {
    if (q_delete.size() == max_deletions || (q_delete.size() >= min_deletions && contrib >= opts.contrib_thresh)) {
      break;
    }

    if (opts.contrib_thresh != 0 && contrib >= opts.contrib_thresh) {
      logger.warn("deleted Q vector with max. contribution = ", contrib);
    } else {
      logger.debug("deleted Q vector with max. contribution = ", contrib);
    }

    q_delete.emplace_back(idx);
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
                      const std::vector<std::size_t>& q_delete, const value_type_abs norm_thresh,
                      const value_type_abs svd_thresh, array::ArrayHandler<Q, Q>& handler, Logger& logger) {
  const auto dims = xspace.dimensions();
  const auto overlap = xspace.data.at(subspace::EqnData::S);
  auto solutions_proj = dspace::construct_projected_solution(solutions, dims, q_delete, logger);
  auto overlap_proj = dspace::construct_projected_solutions_overlap(solutions_proj, overlap, dims, q_delete, logger);
  dspace::remove_null_norm_and_normalise(solutions_proj, overlap_proj, norm_thresh, logger);
  solutions_proj = dspace::remove_null_projected_solutions(solutions_proj, overlap_proj, svd_thresh, logger);
  overlap_proj = dspace::construct_projected_solutions_overlap(solutions_proj, overlap, dims, q_delete, logger);
  dspace::remove_null_norm_and_normalise(solutions_proj, overlap_proj, norm_thresh, logger);
#ifndef NDEBUG
  // a decomposition of the whole subspace overlap, wanted only by this check
  auto overlap_full_subspace = dspace::construct_full_subspace_overlap(solutions_proj, dims, q_delete, overlap, 0);
  auto svd_vecs = svd_system(overlap_full_subspace.rows(), overlap_full_subspace.cols(),
                             array::Span(&overlap_full_subspace(0, 0), overlap_full_subspace.size()), svd_thresh, true);
  assert(svd_vecs.empty() && "P+Q+D subspace should be stable by construction");
#endif
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
    if (norm < norm_thresh) {
      logger.warn(
          std::format("construct_dspace: skipping normalisation of D vector {} with near-zero norm = {:.2e}", i, norm));
      continue;
    }
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
 * @param rx_overlap overlaps <r_j|x_i> of the R space parameters (rows) with the P+Q+D parameters (columns), as
 * already calculated by append_overlap_with_r()
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
                           const subspace::Matrix<value_type>& rx_overlap, const subspace::Dimensions& dims,
                           const CVecRef<P>& pparams, const CVecRef<Q>& qparams, const CVecRef<Q>& dparams,
                           const value_type_abs norm_thresh,
                           ArrayHandlers<R, Q, P>& handlers, Logger& logger) {
  logger.trace("modified_gram_schmidt()");
  const auto nR = rparams.size(), nP = pparams.size(), nQ = qparams.size(), nD = dparams.size();
  // std::cout << "modified_gram_schmidt() nR="<<nR<<" nQ="<<nQ<<" nD="<<nD<<std::endl;
  assert(nP == dims.nP && nQ == dims.nQ && nD == dims.nD);
  assert(rx_overlap.rows() == nR && rx_overlap.cols() == dims.nX());
  // Project the new parameters out of the span of P+Q+D. The existing vectors need not be orthogonal to each other
  // (for example user-supplied initial guesses), so subtract V S^{-1} V^dagger r rather than projecting on each vector
  // separately.
  auto prof = molpro::Profiler::single();
  prof->start("orthoganalise");
  const auto nX = nP + nQ + nD;
  if (nR > 0 and nX > 0) {
    using EMatrix = Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>;
    EMatrix S(nX, nX);
    for (size_t i = 0; i < nX; ++i)
      for (size_t j = 0; j < nX; ++j)
        S(i, j) = overlap(i, j);
    // equilibrate, so that vectors of very different lengths do not make S needlessly ill-conditioned
    Eigen::Vector<value_type, Eigen::Dynamic> scale(nX);
    for (size_t i = 0; i < nX; ++i)
      scale(i) = std::abs(S(i, i)) > 0 ? 1 / std::sqrt(std::abs(S(i, i))) : 0;
    const auto Sinv = (scale.asDiagonal() * S * scale.asDiagonal()).eval().completeOrthogonalDecomposition();
    // the projection needs <x_i|r_j>, the conjugate of the overlaps already calculated
    EMatrix dots(nX, nR);
    for (size_t i = 0; i < nX; ++i)
      for (size_t j = 0; j < nR; ++j)
        dots(i, j) = molpro::linalg::conjugate(rx_overlap(j, i));
    const EMatrix y = scale.asDiagonal() * Sinv.solve(scale.asDiagonal() * dots);
    auto subtract = [&y, &rparams, nR](const auto& xparams, auto& handler, size_t o, size_t n) {
      if (n == 0)
        return;
      subspace::Matrix<value_type> coeff({n, nR});
      for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < nR; ++j)
          coeff(i, j) = -y(o + i, j);
      handler.gemm_outer(coeff, xparams, rparams);
    };
    subtract(pparams, handlers.rp(), dims.oP(), nP);
    subtract(qparams, handlers.rq(), dims.oQ(), nQ);
    subtract(dparams, handlers.rq(), dims.oD(), nD);
  }
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
template <class R, class Q, class P>
auto propose_rspace(IterativeSolver<R, Q, P>& solver, const VecRef<R>& parameters, const VecRef<R>& residuals,
                    subspace::IXSpace<R, Q, P>& xspace, subspace::ISubspaceSolver<R, Q, P>& subspace_solver,
                    ArrayHandlers<R, Q, P>& handlers, Logger& logger,
                    const RSpaceOptions<typename array::ArrayHandler<R, R>::value_type_abs>& r_opts,
                    const QSpaceOptions& q_opts, molpro::profiler::Profiler& profiler) {
  // auto prof = profiler.push("itsolv::propose_rspace"); // FIXME two separate profilers
  auto prof = molpro::Profiler::single();
  logger.trace("itsolv::detail::propose_rspace");
  logger.trace("dimensions {nP, nQ, nD, nW} = ", xspace.dimensions().nP, xspace.dimensions().nQ, xspace.dimensions().nD,
               solver.working_set().size());
  profiler.start("itsolv::ISubspaceSolver::solutions");
  auto solutions = subspace_solver.solutions();
  profiler.stop();
  auto q_delete = limit_qspace_size(xspace.dimensions(), q_opts, solutions, logger);
  logger.debug("delete Q parameter indices = ", q_delete);
  if (!q_delete.empty()) {
    auto prof = profiler.push("construct_dspace");
    auto [dparams, dactions] =
        construct_dspace(solutions, xspace, q_delete, r_opts.norm_thresh, r_opts.svd_thresh, handlers.qq(), logger);
    std::sort(begin(q_delete), end(q_delete), std::greater<int>());
    for (auto iq : q_delete)
      xspace.eraseq(iq);
    auto wdparams = wrap(dparams);
    auto wdactions = wrap(dactions);
    xspace.update_dspace(wdparams, wdactions);
    auto eigenvalues_ref = subspace_solver.eigenvalues();
    subspace_solver.solve(xspace, solutions.rows());
    const auto& eigenvalues_new = subspace_solver.eigenvalues();
    const auto n_compare = std::min(eigenvalues_ref.size(), eigenvalues_new.size());
    auto eigval_error = std::vector<double>{};
    std::transform(std::begin(eigenvalues_ref), std::begin(eigenvalues_ref) + n_compare, std::begin(eigenvalues_new),
                   std::back_inserter(eigval_error), [](auto& e_ref, auto& e_new) { return std::abs(e_ref - e_new); });
    logger.debug("eigenvalue error due to new D space = ", eigval_error);
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
      redundant_parameters(full_overlap, xspace.dimensions().nX(), wresidual.size(), r_opts.svd_thresh, logger);
  prof->stop();
  logger.debug("redundant indices = ", redundant_indices);
  // overlaps of the remaining R parameters with P+Q+D, to save modified_gram_schmidt() from recalculating them
  const auto nX = xspace.dimensions().nX();
  auto kept = std::vector<size_t>{};
  for (size_t j = 0; j < wresidual.size(); ++j)
    if (std::find(redundant_indices.begin(), redundant_indices.end(), int(j)) == redundant_indices.end())
      kept.push_back(j);
  auto rx_overlap = subspace::Matrix<typename decltype(full_overlap)::value_type>({kept.size(), nX});
  for (size_t j = 0; j < kept.size(); ++j)
    for (size_t i = 0; i < nX; ++i)
      rx_overlap(j, i) = full_overlap(nX + kept[j], i);
  util::delete_parameters(redundant_indices, wresidual);
  profiler.start("modified_gram_schmidt");
  prof->start("modified_gram_schmidt");
  auto null_param_indices =
      modified_gram_schmidt(wresidual, xspace.data.at(subspace::EqnData::S), rx_overlap, xspace.dimensions(),
                            xspace.cparamsp(), xspace.cparamsq(), xspace.cparamsd(), r_opts.norm_thresh, handlers,
                            logger);
  profiler.stop();
  prof->stop();
  // Now that there is SVD null_param_indices should always be empty
  logger.debug("null parameters = ", null_param_indices);
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
