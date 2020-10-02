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

//! Proposes an orthonormal set of vectors, removing any that are linearly dependent
//! @returns vector order and linear transformation in that order
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
  return std::tuple<decltype(lin_trans), decltype(norm)>{lin_trans, norm};
}

/*!
 * @brief Construct an orthonormal set from params and a linear transformation matrix
 * @param params vectors to orthonormalise
 * @param lin_trans linear transformation to orthogonal set. Assumed to be from Gram-Schmidt orthogonalisation in lower
 * diagonal form with 1's on diagonal
 * @param norm estimated norm of orthogonalised vectors
 * @param norm_thresh threshold for normalisation at the end
 */
template <class R, typename value_type, typename value_type_abs>
void construct_orthonormal_set(VecRef<R>& params, const subspace::Matrix<value_type>& lin_trans,
                               const std::vector<value_type_abs>& norm, array::ArrayHandler<R, R>& handler,
                               const double norm_thresh) {
  // apply transformation
  // normalise
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
size_t propose_rspace(LinearEigensystem<typename QS::R, typename QS::Q, typename QS::P>& solver,
                      std::vector<typename QS::R>& parameters, std::vector<typename QS::R>& residuals, PS& pspace,
                      QS& qspace, RS& rspace, CS& cspace,
                      subspace::XSpace<RS, QS, PS, CS, typename QS::value_type>& xspace,
                      ArrayHandlers<typename QS::R, typename QS::Q, typename QS::P>& handlers, Logger& logger,
                      typename QS::value_type_abs res_norm_thresh = 1.0e-14) {
  using value_type_abs = typename QS::value_type_abs;
  using value_type = typename QS::value_type;
  using R = typename QS::R;
  auto nw = solver.working_set().size();
  auto wres = wrap<R>(residuals.begin(), std::next(residuals.begin(), nw));
  normalise(wres, res_norm_thresh, handlers.rr(), logger);
  auto lin_trans_and_norm =
      propose_orthonormal_set<R, value_type, value_type_abs>(wres, res_norm_thresh, handlers.rr());
  construct_orthonormal_set(wres, std::get<0>(lin_trans_and_norm), std::get<1>(lin_trans_and_norm), handlers.rr(),
                            res_norm_thresh);
  // Construct overlap matrix for P+Q+C+R subspace
  // P+Q+C can be obtained from xspace, overlap with R has to be calcualted
  // Apply GS procedure on the overlap matrix, suggesting removal of Q space vectors
  // Remove suggest Q vectors
  // Apply linear transformation and normalisation
}
} // namespace detail
} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_PROPOSE_RSPACE_H
