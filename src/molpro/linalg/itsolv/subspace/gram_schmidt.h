#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_GRAM_SCHMIDT_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_GRAM_SCHMIDT_H
#include <molpro/linalg/array/ArrayHandler.h>
#include <molpro/linalg/itsolv/subspace/Matrix.h>
#include <molpro/linalg/itsolv/wrap.h>

namespace molpro::linalg::itsolv::subspace::util {
/*!
 * @brief Performs Gram-Schmidt orthogonalisation without normalisation
 *
 * Any orthogonal vector with a norm less than threshold is not orthogonalised against.
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
 * @param norm_thresh generated vector with norm less than threshold is not used for subsequent orthogonalisation
 * @returns norms of transformed vectors
 */
template <typename T>
std::vector<T> gram_schmidt(const Matrix<T>& s, Matrix<T>& l, double norm_thresh = 1.0e-14) {
  assert(s.rows() == s.cols());
  auto n = s.rows();
  l.fill(0);
  l.resize({n, n});
  auto norm = std::vector<double>(n, 0);
  auto w = std::vector<double>{};
  for (size_t i = 0; i < n; ++i) {
    w.assign(i, 0.);
    for (size_t j = 0; j < i; ++j) {
      for (size_t k = 0; k <= j; ++k) {
        w[j] += s(i, k) * l(j, k);
      }
    }
    for (size_t j = 0; j < i; ++j) {
      if (norm[j] > norm_thresh) {
        for (size_t k = 0; k <= j; ++k) {
          l(i, k) -= w[j] / norm[j] * l(j, k);
        }
      }
    }
    l(i, i) = 1.;
    for (size_t j = 0; j <= i; ++j) {
      for (size_t k = 0; k < j; ++k) {
        norm[i] += 2 * l(i, j) * l(i, k) * s(j, k);
      }
      norm[i] += l(i, j) * l(i, j) * s(j, j);
    }
  }
  std::transform(begin(norm), end(norm), begin(norm), [](auto el) { return std::sqrt(std::abs(el)); });
  return norm;
}

// FIXME This is implicitly constructed in Gram Schmidt and we can make it an opitonal return variable
/*!
 * @brief Construct Gram-Schmidt linear transformation in orthogonal vectors
 *
 * Let {v} be the original vectors, and {u} their orthogonal set.
 *
 * We can construct {u} from {v} using lower triangular linear transformation matrix L
 * u_0 = v_0
 * u_1 = v_1 + L_{1,0} v_0
 * u_i = \sum_{j=0}^i L_{i,j} v_j
 *
 * However, we can also construct u_i using v_i and previous {u}, using a lower triangular transformation matrix T,
 * u_i = v_i - \sum_{j=0}^{i-1} <v_i,u_j> / <u_j, u_j> u_j
 *     = v_i + \sum_j=0^{i-1} T_{i,j} u_j
 *
 * T_{i,j} = - <v_i,u_j> / <u_j, u_j>
 *         = - \sum_{k=0}^{j} L_{j,k} <v_i,v_k> / <u_j, u_j>
 *
 * T_{i,j} = 0 if i >= j
 *
 * @param overlap Overlap matrix of original non-orthogonal vectors
 * @param lin_trans Gram-Schmidt linear transformation to an orthogonal set in terms of original vectors
 * @param norm norms of orthogonal vectors constructed from lin_trans
 * @Returns linear transformation set
 * @return
 */
template <typename value_type, typename value_type_abs>
Matrix<value_type> construct_lin_trans_in_orthogonal_set(const Matrix<value_type>& overlap,
                                                         const Matrix<value_type>& lin_trans,
                                                         const std::vector<value_type_abs>& norm) {
  const auto nrows = lin_trans.rows();
  const auto ncols = lin_trans.cols();
  auto t = Matrix<value_type>({nrows, ncols});
  for (size_t i = 0; i < nrows; ++i) {
    for (size_t j = 0; j < i; ++j) {
      for (size_t k = 0; k <= j; ++k) {
        t(i, j) -= lin_trans(j, k) * overlap(i, k) / std::pow(norm[j], 2);
      }
    }
  }
  return t;
}

/*!
 * @brief Apply modified Gram-Schmidt procedure to orthonormalise parameters
 *
 * New vectors with norm less than threshold are considered null and are not normalised.
 * Their indices in params are returned.
 *
 * @tparam R
 * @tparam value_type_abs
 * @param params
 * @param handler
 * @param null_thresh
 * @return indices of parameters that became null after GS procedure
 */
template <class R, typename value_type_abs>
auto modified_gram_schmidt(VecRef<R>& params, array::ArrayHandler<R, R>& handler, value_type_abs null_thresh) {
  auto null_param_indices = std::vector<size_t>{};
  const size_t n = params.size();
  for (size_t i = 0; i < n; ++i) {
    auto norm = handler.dot(params[i], params[i]);
    norm = std::sqrt(std::abs(norm));
    if (norm > null_thresh) {
      handler.scal(1. / norm, params[i]);
      for (size_t j = i + 1; j < n; ++j) {
        auto ov = handler.dot(params[i], params[j]);
        handler.axpy(-ov, params[i], params[j]);
      }
    } else {
      null_param_indices.emplace_back(i);
    }
  }
  return null_param_indices;
}

} // namespace molpro::linalg::itsolv::subspace::util

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_GRAM_SCHMIDT_H
