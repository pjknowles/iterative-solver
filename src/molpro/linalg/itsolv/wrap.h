#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_WRAP_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_WRAP_H
#include <functional>
#include <limits>
#include <type_traits>
#include <vector>

namespace molpro {
namespace linalg {
namespace itsolv {

template <class A>
using VecRef = std::vector<std::reference_wrapper<A>>;

template <class A>
using CVecRef = std::vector<std::reference_wrapper<const A>>;

template <class T>
struct decay {
  using type = std::decay_t<T>;
};

template <class T>
struct decay<std::reference_wrapper<T>> {
  using type = std::decay_t<T>;
};

template <class T>
using decay_t = typename decay<T>::type;

//! Takes a vector of containers and returns a vector of references to each element
template <class R>
auto wrap(const std::vector<R>& vec) {
  auto w = CVecRef<decay_t<R>>{};
  std::copy(begin(vec), end(vec), std::back_inserter(w));
  return w;
}

//! Takes a vector of containers and returns a vector of references to each element
template <class R>
auto wrap(std::vector<R>& vec) {
  auto w = VecRef<decay_t<R>>{};
  std::copy(begin(vec), end(vec), std::back_inserter(w));
  return w;
}

//! Takes a map of containers and returns a vector of references to each element in the same order
template <class R>
auto wrap(const std::map<size_t, R>& vec) {
  auto w = CVecRef<decay_t<R>>{};
  std::transform(begin(vec), end(vec), std::back_inserter(w), [](const auto& v) { return v.second; });
  return w;
}

//! Takes a map of containers and returns a vector of references to each element in the same order
template <class R>
auto wrap(std::map<size_t, R>& vec) {
  auto w = VecRef<decay_t<R>>{};
  std::transform(begin(vec), end(vec), std::back_inserter(w), [](const auto& v) { return v.second; });
  return w;
}

} // namespace itsolv
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_WRAP_H
