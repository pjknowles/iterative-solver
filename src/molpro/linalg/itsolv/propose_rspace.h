#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_PROPOSE_RSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_PROPOSE_RSPACE_H
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/subspace/QSpace.h>
#include <molpro/linalg/itsolv/subspace/XSpaceI.h>
#include <molpro/linalg/itsolv/subspace/gram_schmidt.h>
#include <molpro/linalg/itsolv/subspace/util.h>

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
auto propose_orthonormal_set(VecRef<R>& params, const double norm_thresh, array::ArrayHandler<R, R>& handler,
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
  return std::tuple<decltype(lin_trans), decltype(norm)>{lin_trans, norm};
}
/*!
 * @brief Uses overlap matrix to construct an orthogonal set of R params, removing any Q parameters that are linearly
 *   dependent
 * @param overlap overlap matrix of current subspace + initial R parameters
 * @param qspace qspace container
 * @param oQ offset to the Q parameters in the full subspace
 * @param nW number of working parameters
 * @param res_norm_thresh norm threshold for Gram-Schmidt procedure
 * @return
 */
template <class R, class Q, class P, typename value_type, typename value_type_abs>
auto prepare_orthogonal_rset(subspace::Matrix<value_type>& overlap, subspace::XSpaceI<R, Q, P>& xspace, const size_t nW,
                             Logger& logger, value_type_abs res_norm_thresh) {
  logger.msg("prepare_orthogonal_rset()", Logger::Trace);
  auto norm = std::vector<value_type_abs>{};
  auto lin_trans = subspace::Matrix<value_type>{};
  for (bool done = false; !done;) {
    const auto nQ = xspace.dimensions().nQ;
    const auto oQ = xspace.dimensions().oQ;
    const auto oN = oQ + nQ;
    norm = subspace::util::gram_schmidt(overlap, lin_trans, res_norm_thresh);
    auto it = std::find_if(std::next(begin(norm), oN), end(norm),
                           [res_norm_thresh](auto el) { return el < res_norm_thresh; });
    done = (it == end(norm) || nQ == 0);
    if (!done) {
      auto i = std::distance(begin(norm), it);
      logger.msg("parameter index i = " + std::to_string(i) + " norm = " + std::to_string(*it), Logger::Info);
      auto normalised_overlap = std::vector<value_type>{};
      for (size_t j = 0; j < nQ; ++j) {
        normalised_overlap.emplace_back(std::abs(overlap(i, oQ + j)) / std::sqrt(std::abs(overlap(oQ + j, oQ + j))));
      }
      auto it_max = std::max_element(begin(normalised_overlap), end(normalised_overlap));
      auto iq_erase = std::distance(begin(normalised_overlap), it_max);
      const auto ix = oQ + iq_erase;
      logger.msg("removing parameter index = " + std::to_string(ix), Logger::Info);
      xspace.eraseq(iq_erase);
      overlap.remove_row_col(ix, ix);
    }
  }
  return std::tuple<decltype(lin_trans), decltype(norm)>{lin_trans, norm};
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
 * @param oP offset to the start of P parameters
 * @param oQ offset to the start of Q parameters
 * @param oC offset to the start of C parameters
 * @param oN offset to the start of new parameters
 * @param handlers array handlers
 * @param logger logger
 */
template <class R, class Q, class P, typename value_type>
auto append_overlap_with_r(const subspace::Matrix<value_type>& overlap, const CVecRef<R>& params,
                           const CVecRef<P>& pparams, const CVecRef<Q>& qparams,
                           const subspace::xspace::Dimensions& dims, ArrayHandlers<R, Q, P>& handlers, Logger& logger) {
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
 * @param cparams C space parameters
 * @param oP offset to start of P space paramters
 * @param oQ offset to start of Q space paramters
 * @param oC offset to start of C space paramters
 * @param oN offset to start of the new parameters
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

/*!
 * \brief Proposes new parameters for the subspace from the preconditioned residuals.
 *
 * Outline
 * -------
 * Basic procedure:
 *  - Gram-schmidt orthonormalise residuals amongst themselves
 *  - Gram-schmidt orthogonalise residuals against the old Q space and current solutions
 *  - Ensure that the resultant Q space is not linearly dependent
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
std::vector<unsigned int> propose_rspace(LinearEigensystem<R, Q, P>& solver, std::vector<R>& parameters,
                                         std::vector<R>& residuals, subspace::XSpaceI<R, Q, P>& xspace,
                                         ArrayHandlers<R, Q, P>& handlers, Logger& logger,
                                         value_type_abs res_norm_thresh = 1.0e-6) {
  using value_type = typename std::decay_t<decltype(xspace)>::value_type;
  logger.msg("itsolv::detail::propose_rspace", Logger::Trace);
  auto nW = solver.working_set().size();
  auto wresidual = wrap<R>(residuals.begin(), std::next(residuals.begin(), nW));
  normalise(wresidual, handlers.rr(), logger);
  auto [lin_trans, norm] =
      propose_orthonormal_set<R, value_type, value_type_abs>(wresidual, res_norm_thresh, handlers.rr(), logger);
  nW = wresidual.size();
  auto new_indices = find_ref(wresidual, begin(residuals), end(residuals));
  auto new_working_set = std::vector<unsigned int>{};
  auto wparams = wrap<R>(parameters.begin(), parameters.begin() + nW);
  for (auto i : new_indices) {
    new_working_set.emplace_back(solver.working_set().at(i));
  }
  construct_orthonormal_set(wresidual, lin_trans, norm, handlers.rr());
  auto ov = append_overlap_with_r(xspace.data.at(subspace::EqnData::S), cwrap(wresidual), xspace.cparamsp(),
                                  xspace.cparamsq(), xspace.dimensions(), handlers, logger);
  std::tie(lin_trans, norm) = prepare_orthogonal_rset(ov, xspace, nW, logger, res_norm_thresh);
  if (logger.data_dump) {
    logger.msg("full overlap = " + subspace::as_string(ov), Logger::Info);
    logger.msg("linear transformation = " + subspace::as_string(lin_trans), Logger::Info);
    logger.msg("norm = ", norm.begin(), norm.end(), Logger::Info);
  }
  construct_orthonormal_Rparams(wparams, wresidual, lin_trans, norm, xspace.cparamsp(), xspace.cparamsq(), handlers);
  normalise(wparams, handlers.rr(), logger);
  return new_working_set;
}
} // namespace molpro::linalg::itsolv::detail

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_PROPOSE_RSPACE_H
