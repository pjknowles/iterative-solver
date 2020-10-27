#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_RESET_DSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_RESET_DSPACE_H
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/subspace/QSpace.h>
#include <molpro/linalg/itsolv/subspace/XSpaceI.h>
#include <molpro/linalg/itsolv/subspace/gram_schmidt.h>
#include <molpro/linalg/itsolv/subspace/util.h>

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
 * @param parameters parameters to propose as R space
 * @param xspace current X space container
 * @param solutions current solutions in the subspace stored as rows
 * @param norm_thresh threshold for selecting whether the full solution can be added as part of resetting
 * @param handlers array handlers
 * @param logger logger
 * @return new working set, only indicating number of parameters that need their action evaluated
 */
template <class R, class Q, class P, typename value_type, typename value_type_abs>
auto reset_dspace(LinearEigensystem<R, Q, P>& solver, std::vector<R>& parameters, subspace::XSpaceI<R, Q, P>& xspace,
                  const subspace::Matrix<value_type>& solutions, value_type_abs norm_thresh,
                  ArrayHandlers<R, Q, P>& handlers, Logger& logger) {
  auto overlap = construct_overlap_with_solutions(solutions, xspace.data.at(subspace::EqnData::S), xspace.dimensions());
  auto lin_trans = subspace::Matrix<value_type>{};
  auto norm = subspace::util::gram_schmidt(overlap, lin_trans, 0);
  const auto nC = solutions.rows();
  const auto nP = xspace.dimensions().nP;
  const auto nQ = xspace.dimensions().nQ;
  const auto nPQ = nP + nQ;
  auto lin_trans_orth_sol = subspace::Matrix<value_type>({nC, nPQ + nC});
  lin_trans_orth_sol.slice() = lin_trans.slice({nPQ, 0}, lin_trans.dimensions());
  auto norm_orth_sol = std::vector<value_type_abs>(std::begin(norm) + nPQ, std::begin(norm) + nPQ + nC);
  lin_trans_orth_sol = transform_PQSol_to_PQD_subspace(lin_trans_orth_sol, solutions, xspace.dimensions());
  const auto nR = std::min(parameters.size(), xspace.dimensions().nD);
  const auto nDnew = nC - nR;
  subspace::Matrix<value_type> lin_trans_R, lin_trans_D;
  std::tie(lin_trans_R, lin_trans_D) =
      propose_R_and_D_params(lin_trans_orth_sol, norm_orth_sol, solutions, nR, nDnew, norm_thresh);
  std::vector<Q> dparams, dactions;
  std::tie(dparams, dactions) =
      construct_R_and_D_params(parameters, lin_trans_R, lin_trans_D, xspace.cparamsp(), xspace.cactionsp(),
                               xspace.cparamsq(), xspace.cactionsq(), xspace.cparamsd(), xspace.cactionsd(), handlers);
  auto wdparams = wrap(dparams);
  auto wdactions = wrap(dactions);
  xspace.update_dspace(wdparams, wdactions, subspace::Matrix<value_type>({nDnew, nR}));
  const auto& working_set = solver.working_set();
  auto new_working_set = std::vector<unsigned int>(std::begin(working_set), std::begin(working_set) + nR);
  return new_working_set;
}

} // namespace detail
} // namespace itsolv
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_RESET_DSPACE_H
