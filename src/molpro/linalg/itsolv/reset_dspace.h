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
      max_size_qspace += nD;
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
 * @brief Removes part of D space and uses it as R space so that the exact action can be calculated
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
  // Construct overlap matrix P+Q+C
  auto ov = construct_overlap_with_solutions(solutions, xspace.data.at(subspace::EqnData::S), xspace.dimensions());
  // Orthogonalise C against P+Q
  // Construct new R and left-over D parameters
  auto new_working_set = solver.working_set();
  return new_working_set;
}
} // namespace detail
} // namespace itsolv
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_RESET_DSPACE_H
