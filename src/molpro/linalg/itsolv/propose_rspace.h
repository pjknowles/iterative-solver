#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_PROPOSE_RSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_PROPOSE_RSPACE_H
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/subspace/QSpace.h>
#include <molpro/linalg/itsolv/subspace/gram_schmidt.h>
#include <molpro/linalg/itsolv/subspace/util.h>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace detail {

template <class R>
void normalise(VecRef<R>& params, double thresh, array::ArrayHandler<R, R>& handler, Logger& logger) {
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
auto propose_orthonormal_set(VecRef<R>& params, const double norm_thresh, array::ArrayHandler<R, R>& handler) {
  auto lin_trans = subspace::Matrix<value_type>{};
  auto norm = std::vector<value_type_abs>{};
  auto ov = subspace::util::overlap(params, handler);
  for (auto done = false; !done;) {
    norm = subspace::util::gram_schmidt(ov, lin_trans);
    auto it_small = std::find_if(begin(norm), end(norm), [norm_thresh](const auto& el) { return el < norm_thresh; });
    done = (it_small == end(norm));
    if (!done) {
      auto i = std::distance(begin(norm), it_small);
      auto it_remove = std::next(begin(params), i);
      params.erase(it_remove);
      ov.remove_row_col(i, i);
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
template <class QS, typename value_type, typename value_type_abs>
auto prepare_orthogonal_rset(subspace::Matrix<value_type>& overlap, QS& qspace, const size_t oQ, const size_t nW,
                             value_type_abs res_norm_thresh) {
  auto norm = std::vector<value_type_abs>{};
  auto lin_trans = subspace::Matrix<value_type>{};
  for (bool done = false; !done;) {
    const auto nX = overlap.rows() - nW;
    const auto nQ = qspace.size();
    norm = subspace::util::gram_schmidt(overlap, lin_trans);
    auto it = std::find_if(std::next(begin(norm), nX), end(norm),
                           [res_norm_thresh](auto el) { return el < res_norm_thresh; });
    done = (it == end(norm) || nQ == 0);
    if (!done) {
      auto i = std::distance(begin(norm), it);
      // FIXME Not sure whether to just use the overlap or linear transformation matrix
      auto normalised_transform = std::vector<value_type>{};
      for (size_t j = 0; j < nQ; ++j) {
        normalised_transform.emplace_back(lin_trans(i, oQ + j) / std::sqrt(std::abs(overlap(oQ + j, oQ + j))));
      }
      auto it_max = std::max_element(begin(normalised_transform), end(normalised_transform));
      auto iq_erase = std::distance(begin(normalised_transform), it_max);
      qspace.erase(iq_erase);
      overlap.remove_row_col(oQ + iq_erase, oQ + iq_erase);
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
 * @param ov overlap matrix
 * @param params new parameters
 * @param pparams P space parameters
 * @param qparams Q space parameters
 * @param cparams C space parameters
 * @param oP offset to the start of P parameters
 * @param oQ offset to the start of Q parameters
 * @param oC offset to the start of C parameters
 * @param oN offset to the start of new parameters
 * @param handlers array handlers
 * @param logger logger
 */
template <class R, class Q, class P, typename value_type>
auto append_overlap_with_r(subspace::Matrix<value_type> ov, VecRef<R>& params, const CVecRef<P>& pparams,
                           const CVecRef<Q>& qparams, const CVecRef<Q>& cparams, const size_t oP, const size_t oQ,
                           const size_t oC, const size_t oN, ArrayHandlers<R, Q, P>& handlers, Logger& logger) {
  auto nP = pparams.size();
  auto nQ = qparams.size();
  auto nC = cparams.size();
  auto nN = params.size();
  ov.resize({ov.rows() + nN, ov.cols() + nN});
  ov.slice({oN, oN}, {oN + nN, oN + nN}) = subspace::util::overlap(params, handlers.rr());
  ov.slice({oN, oP}, {oN + nN, oP + nP}) = subspace::util::overlap(params, pparams, handlers.rp());
  ov.slice({oN, oQ}, {oN + nN, oQ + nQ}) = subspace::util::overlap(params, qparams, handlers.rq());
  ov.slice({oN, oC}, {oN + nN, oC + nC}) = subspace::util::overlap(params, cparams, handlers.rq());
  auto copy_upper_to_lower = [&ov, oN, nN](size_t oX, size_t nX) {
    for (size_t i = 0; i < nX; ++i)
      for (size_t j = 0; j < nN; ++j)
        ov(oX + i, oN + j) = ov(oN + j, oX + i);
  };
  copy_upper_to_lower(oP, nP);
  copy_upper_to_lower(oQ, nQ);
  copy_upper_to_lower(oC, nC);
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
                                   const CVecRef<Q>& qparams, const CVecRef<Q>& cparams, const size_t oP,
                                   const size_t oQ, const size_t oC, const size_t oN,
                                   ArrayHandlers<R, Q, P>& handlers) {
  assert(params.size() == residuals.size());
  for (size_t i = 0; i < params.size(); ++i) {
    handlers.rr().copy(params.at(i), residuals.at(i));
  }
  for (size_t i = 0; i < params.size(); ++i) {
    for (size_t j = 0; j < pparams.size(); ++j) {
      handlers.rp().axpy(lin_trans(oN + i, oP + j), pparams.at(j), params.at(i));
    }
    for (size_t j = 0; j < qparams.size(); ++j) {
      handlers.rq().axpy(lin_trans(oN + i, oQ + j), qparams.at(j), params.at(i));
    }
    for (size_t j = 0; j < cparams.size(); ++j) {
      handlers.rq().axpy(lin_trans(oN + i, oC + j), cparams.at(j), params.at(i));
    }
    for (size_t j = 0; j < i; ++j) {
      handlers.rr().axpy(lin_trans(oN + i, oN + j), residuals.at(j), params.at(i));
    }
  }
  for (size_t i = 0; i < params.size(); ++i) {
    handlers.rr().scal(1. / norm.at(i), params.at(i));
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
template <class PS, class QS, class RS, class CS>
std::vector<unsigned int> propose_rspace(LinearEigensystem<typename QS::R, typename QS::Q, typename QS::P>& solver,
                                         std::vector<typename QS::R>& parameters,
                                         std::vector<typename QS::R>& residuals, PS& pspace, QS& qspace, RS& rspace,
                                         CS& cspace, subspace::XSpace<RS, QS, PS, CS, typename QS::value_type>& xspace,
                                         ArrayHandlers<typename QS::R, typename QS::Q, typename QS::P>& handlers,
                                         Logger& logger, typename QS::value_type_abs res_norm_thresh = 1.0e-14) {
  using value_type_abs = typename QS::value_type_abs;
  using value_type = typename QS::value_type;
  using R = typename QS::R;
  logger.msg("itsolv::detail::propose_rspace", Logger::Trace);
  auto nW = solver.working_set().size();
  auto wresidual = wrap<R>(residuals.begin(), std::next(residuals.begin(), nW));
  normalise(wresidual, res_norm_thresh, handlers.rr(), logger);
  auto lin_trans = subspace::Matrix<value_type>{};
  auto norm = std::vector<value_type_abs>{};
  std::tie(lin_trans, norm) =
      propose_orthonormal_set<R, value_type, value_type_abs>(wresidual, res_norm_thresh, handlers.rr());
  nW = wresidual.size();
  auto new_indices = find_ref(wresidual, begin(residuals), end(residuals));
  auto new_working_set = std::vector<unsigned int>{};
  auto wparams = wrap<R>(parameters.begin(), parameters.begin() + nW);
  for (auto i : new_indices) {
    new_working_set.emplace_back(solver.working_set().at(i));
  }
  construct_orthonormal_set(wresidual, lin_trans, norm, handlers.rr());
  xspace.build_subspace(rspace, qspace, pspace, cspace);
  auto ov = append_overlap_with_r(xspace.data.at(subspace::EqnData::S), wresidual, pspace.cparams(), qspace.cparams(),
                                  cspace.cparams(), xspace.dimensions().oP, xspace.dimensions().oQ,
                                  xspace.dimensions().oC, xspace.dimensions().nX, handlers, logger);
  std::tie(lin_trans, norm) = prepare_orthogonal_rset(ov, qspace, xspace.dimensions().oQ, nW, res_norm_thresh);
  construct_orthonormal_Rparams(wparams, wresidual, lin_trans, norm, pspace.cparams(), qspace.cparams(),
                                cspace.cparams(), xspace.dimensions().oP, xspace.dimensions().oQ,
                                xspace.dimensions().oC, xspace.dimensions().nX, handlers);
  normalise(wparams, res_norm_thresh, handlers.rr(), logger);
  return new_working_set;
}
} // namespace detail
} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_PROPOSE_RSPACE_H
