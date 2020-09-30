#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_CHECK_CONDITIONING_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_CHECK_CONDITIONING_H
#include <molpro/linalg/itsolv/helper.h>
#include <molpro/linalg/itsolv/subspace/PSpace.h>
#include <molpro/linalg/itsolv/subspace/QSpace.h>
#include <molpro/linalg/itsolv/subspace/RSpace.h>
#include <molpro/linalg/itsolv/subspace/XSpace.h>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace subspace {
namespace xspace {
namespace detail {
/*!
 * @brief Performs Gram-Schmidt orthogonalisation
 *
 * For a vector set {v} with overlap matrix \f$ S_{ij}= v_i . v_j \f$, generate a linear transformation L to a new
 * orthogonal vector set \f$ u_i = \sum_j L_{ij} v_j \f$.
 *
 * u_0 = v_0
 * u_1 = v_1 - <v_1, u_0> / <u_0, u_0>
 * u_i = v_i -\sum_j=0^{i-1} <v_i, u_j> / <u_j, u_j> u_j
 *     = v_i -\sum_j=0^{i-1} <v_i, u_j> / <u_j, u_j> \sum_{k=0}^j L_{jk} v_k
 *     = v_i -\sum_j=0^{i-1} (\sum_{k=0}^{j} W_{ij} / N_{jj} L_{jk}) v_k
 *
 * Let the overlap between u and v be, \f$ W_{ij} = <v_i, u_j> \f$
 *
 * W_{ij} = <v_i, u_j> = \sum_{k=0}^{j} S_{ik} L_{jk}
 * N_{ii} = <u_i, u_i>
 *        = \sum_{k} L_{ik} <v_k, u_i>
 *        = \sum_{km} L_{ik} L_{im} S_{km}
 *
 * L_{0,j} = \delta_{i,j}
 * L_{i,k} = \delta_{i,k} - (1 + \delta_{i,k})(\sum_{j=0}^{k} W_{ij} / N_{jj} L_{jk})
 *
 * @param s overlap matrix
 * @param l row matrix with linear transformation into orthonormal basis. Upper triangular component is zero.
 * @returns norms of transformed vectors
 */
template <typename T>
std::vector<T> gram_schmidt(const Matrix<T>& s, Matrix<T>& l) {
  assert(s.rows() == s.cols());
  auto n = s.rows();
  l.fill(0);
  l.resize({n, n});
  auto norm = std::vector<double>{};
  auto w = std::vector<double>{};
  for (size_t i = 0; i < n; ++i) {
    w.assign(i, 0.);
    for (size_t j = 0; j < i; ++j) {
      for (size_t k = 0; k <= j; ++k) {
        w[j] += s(i, k) * l(j, k);
      }
    }
    for (size_t j = 0; j < i; ++j) {
      for (size_t k = 0; k <= j; ++k) {
        l(i, k) -= w[j] / norm[j] * l(j, k);
      }
    }
    l(i, i) = 1.;
    norm.emplace_back(0);
    for (size_t j = 0; j <= i; ++j) {
      for (size_t k = 0; k < j; ++k) {
        norm.back() += 2 * l(i, j) * l(i, k) * s(j, k);
      }
      norm.back() += l(i, j) * l(i, j) * s(j, j);
    }
  }
  std::transform(begin(norm), end(norm), begin(norm), [](auto el) { return std::sqrt(std::abs(el)); });
  return norm;
}

} // namespace detail

/*!
 * @brief Generates linear transformation of X space to an orthogonal subspace and removes any redundant Q vectors.
 * @param xs X space container
 * @param rs R space container
 * @param qs Q space container
 * @param ps P space container
 * @param cs C space container
 * @param lin_trans  Linear transformation matrix to an orthogonal subspace
 * @param norm_threshold
 * @param logger
 */
template <class R, class P, class Q, class ST>
void check_conditioning_gram_schmidt(
    XSpace<RSpace<R, Q, P>, QSpace<R, Q, P>, PSpace<R, P>, CSpace<R, Q, P, ST>, ST>& xs, RSpace<R, Q, P>& rs,
    QSpace<R, Q, P>& qs, PSpace<R, P>& ps, CSpace<R, Q, P, ST> cs, Matrix<ST>& lin_trans, double norm_threshold,
    Logger& logger) {
  logger.msg("xspace::check_conditioning_gram_schmidt", Logger::Trace);
  bool stable = false;
  auto candidates = std::vector<size_t>{qs.size()};
  std::iota(begin(candidates), end(candidates), size_t{xs.dimensions().oQ});
  while (!stable && !candidates.empty()) {
    const auto& dim = xs.dimensions();
    const auto& s = xs.data[EqnData::S];
    auto norm = detail::gram_schmidt(s, lin_trans);
    if (logger.data_dump)
      logger.msg("norm after Gram-Schmidt = ", begin(norm), end(norm), Logger::Info);
    auto imin = std::find_if(begin(candidates), end(candidates),
                             [&norm, norm_threshold](auto c) { return norm[c] < norm_threshold; });
    stable = (imin == candidates.end());
    if (!stable) {
      logger.msg("removing candidate from q space i =" + std::to_string(*imin) +
                     ", norm = " + Logger::scientific(norm[*imin]),
                 Logger::Debug);
      qs.erase(*imin);
      xs.build_subspace(rs, qs, ps, cs);
      candidates.resize(qs.size());
    }
  }
}

} // namespace xspace
} // namespace subspace
} // namespace itsolv
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_CHECK_CONDITIONING_H
