#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_GRAM_SCHMIDT_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_GRAM_SCHMIDT_H
#include <molpro/linalg/itsolv/helper.h>
#include <molpro/linalg/itsolv/subspace/PSpace.h>
#include <molpro/linalg/itsolv/subspace/QSpace.h>
#include <molpro/linalg/itsolv/subspace/RSpace.h>
#include <molpro/linalg/itsolv/subspace/XSpace.h>
namespace molpro {
namespace linalg {
namespace itsolv {
namespace subspace {
namespace util {
/*!
 * @brief Performs Gram-Schmidt orthogonalisation without normalisation
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
        l(i, k) -= w[j] / norm[j] * l(j, k); // FIXME skip if norm is too small, simply set to 0
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
} // namespace util
} // namespace subspace
} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_GRAM_SCHMIDT_H
