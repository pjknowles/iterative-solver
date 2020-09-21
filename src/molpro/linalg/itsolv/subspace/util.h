#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_UTIL_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_UTIL_H
#include <limits>
#include <molpro/linalg/itsolv/ArrayHandlers.h>
#include <molpro/linalg/itsolv/subspace/Matrix.h>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace subspace {
namespace util {

//! Takes a vector of containers and returns a vector of references to each element
template <class R>
auto wrap(const std::vector<R>& vec) {
  auto w = std::vector<std::reference_wrapper<const R>>{};
  std::copy(begin(vec), end(vec), std::back_inserter(w));
  return w;
}

//! Takes a vector of containers and returns a vector of references to each element
template <class R>
auto wrap(std::vector<R>& vec) {
  auto w = std::vector<std::reference_wrapper<R>>{};
  std::copy(begin(vec), end(vec), std::back_inserter(w));
  return w;
}

//! Takes a map of containers and returns a vector of references to each element in the same order
template <class R>
auto wrap(const std::map<size_t, R>& vec) {
  auto w = std::vector<std::reference_wrapper<const R>>{};
  std::transform(begin(vec), end(vec), std::back_inserter(w), [](const auto& v) { return v.second; });
  return w;
}

//! Takes a map of containers and returns a vector of references to each element in the same order
template <class R>
auto wrap(std::map<size_t, R>& vec) {
  auto w = std::vector<std::reference_wrapper<R>>{};
  std::transform(begin(vec), end(vec), std::back_inserter(w), [](const auto& v) { return v.second; });
  return w;
}

//! Takes a map of container references and returns a vector in the same order
template <class R>
std::vector<std::reference_wrapper<R>> wrap(const std::map<size_t, std::reference_wrapper<R>>& vec) {
  auto w = std::vector<std::reference_wrapper<R>>{};
  std::transform(begin(vec), end(vec), std::back_inserter(w), [](const auto& v) { return v.second; });
  return w;
}

//! Calculates overlap matrix between left and right vectors
template <class R, class Q>
Matrix<double> overlap(const std::vector<std::reference_wrapper<R>>& left,
                       const std::vector<std::reference_wrapper<Q>>& right,
                       array::ArrayHandler<std::decay_t<R>, std::decay_t<Q>>& handler) {
  auto m = Matrix<double>({left.size(), right.size()});
  for (size_t i = 0; i < m.rows(); ++i)
    for (size_t j = 0; j < m.cols(); ++j)
      m(i, j) = handler.dot(left[i], right[j]);
  return m;
}

//! Calculates overlap matrix for a parameter set. Matrix is symmetric by construction.
template <class R>
Matrix<double> overlap(const std::vector<std::reference_wrapper<R>>& params,
                       array::ArrayHandler<std::decay_t<R>, std::decay_t<R>>& handler) {
  auto m = Matrix<double>({params.size(), params.size()});
  for (size_t i = 0; i < m.rows(); ++i)
    for (size_t j = 0; j <= i; ++j)
      m(i, j) = m(j, i) = handler.dot(params[i], params[j]);
  return m;
}

template <typename T>
void matrix_symmetrize(Matrix<T>& mat) {
  assert(mat.rows() == mat.cols() && "must be a square matrix");
  for (size_t i = 0; i < mat.rows(); ++i)
    for (size_t j = 0; j < i; ++j)
      mat(i, j) = mat(j, i) = 0.5 * (mat(i, j) + mat(j, i));
}

//! Return maximum element in a matrix along specified rows and columns
template <typename T>
typename Matrix<T>::coord_type max_element_index(const std::list<size_t>& rows, const std::list<size_t>& cols,
                                                 const Matrix<T>& mat) {
  auto max_el = std::numeric_limits<T>::lowest();
  auto ind = typename Matrix<T>::coord_type{0, 0};
  for (auto i : rows) {
    for (auto j : cols) {
      if (mat(i, j) > max_el) {
        max_el = mat(i, j);
        ind = {i, j};
      }
    }
  }
  return ind;
}

/*!
 * @brief Returns order of rows in a matrix slice that brings it closest to identity.
 *
 * @code{.cpp}
 * auto order = eye_order(mat);
 * for( size_t i = 0; i < n_rows; ++i){
 *   mat_new.row(i) = mat.row(order[i]);
 * }
 *
 * @endcode
 * @tparam Slice Slice type
 * @param mat  square matrix
 */
template <typename Slice>
std::vector<size_t> eye_order(const Slice& mat) {
  auto dim = mat.dimensions();
  auto rows = std::list<size_t>{};
  auto cols = std::list<size_t>{};
  for (size_t i = 0; i < dim.first; ++i) {
    rows.emplace_back(i);
    cols.emplace_back(i);
  }
  auto order = std::vector<size_t>(dim.first);
  size_t i, j;
  while (!rows.empty() && !cols.empty()) {
    std::tie(i, j) = max_element_index(rows, cols, mat);
    order.at(j) = i;
    auto it_row = std::find(begin(rows), end(rows), i);
    auto it_col = std::find(begin(cols), end(cols), j);
    rows.erase(it_row);
    cols.erase(it_col);
  }
  return order;
}

} // namespace util
} // namespace subspace
} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_UTIL_H
