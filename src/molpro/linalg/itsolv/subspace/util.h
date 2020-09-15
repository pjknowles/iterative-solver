#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_UTIL_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_UTIL_H
#include <molpro/linalg/itsolv/ArrayHandlers.h>
#include <molpro/linalg/itsolv/subspace/Matrix.h>
namespace molpro {
namespace linalg {
namespace itsolv {
namespace subspace {
namespace util {

//! Takes a vector of containers and returns a vector of references to each element
template <class R>
std::vector<std::reference_wrapper<R>> wrap(const std::vector<R>& vec) {
  auto w = std::vector<std::reference_wrapper<R>>{};
  for (const auto& v : vec) {
    w.emplace_back(v);
  }
  return w;
}
//! Takes a map of containers and returns a vector of references to each element in the same order
template <class R>
std::vector<std::reference_wrapper<R>> wrap(const std::map<size_t, R>& vec) {
  auto w = std::vector<std::reference_wrapper<R>>{};
  for (const auto& v : vec) {
    w.emplace_back(v.second);
  }
  return w;
}

//! Takes a map of container references and returns a vector in the same order
template <class R>
std::vector<std::reference_wrapper<R>> wrap(const std::map<size_t, std::reference_wrapper<R>>& vec) {
  auto w = std::vector<std::reference_wrapper<R>>{};
  for (const auto& v : vec) {
    w.emplace_back(v.second);
  }
  return w;
}

//! Calculates overlap matrix between left and right vectors
template <class R, class Q>
Matrix<double> overlap(const std::vector<std::reference_wrapper<R>>& left,
                       const std::vector<std::reference_wrapper<Q>>& right, const array::ArrayHandler<R, Q>& handler) {
  auto m = Matrix<double>({left.size(), right.size()});
  for (size_t i = 0; i < m.rows(); ++i)
    for (size_t j = 0; j < m.cols(); ++j)
      m(i, j) = handler.dot(left[i], right[j]);
  return m;
}

//! Calculates overlap matrix for a parameter set. Matrix is symmetric by construction.
template <class R>
Matrix<double> overlap(const std::vector<std::reference_wrapper<R>>& params, const array::ArrayHandler<R, R>& handler) {
  auto m = Matrix<double>({params.size(), params.size()});
  for (size_t i = 0; i < m.rows(); ++i)
    for (size_t j = 0; j <= i; ++j)
      m(i, j) = m(j, i) = handler.dot(params[i], params[j]);
  return m;
}

} // namespace util
} // namespace subspace
} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_UTIL_H
