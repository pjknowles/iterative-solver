#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_WRAP_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_WRAP_H
#include <functional>
#include <iterator>
#include <type_traits>
#include <vector>

namespace molpro::linalg::itsolv {

template <class A>
using VecRef = std::vector<std::reference_wrapper<A>>;

template <class A>
using CVecRef = std::vector<std::reference_wrapper<const A>>;

//! Decays CV and reference qualifiers same as std::decay, but also decays std::reference_wrapper to its underlying type
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
  for (const auto& v : vec)
    w.emplace_back(std::cref(v));
  return w;
}

//! Takes a vector of containers and returns a vector of references to each element
template <class R>
auto wrap(std::vector<R>& vec) {
  auto w = VecRef<decay_t<R>>{};
  for (auto& v : vec)
    w.emplace_back(std::ref(v));
  return w;
}

//! Takes a vector of containers and returns a vector of references to each element
template <class R>
auto cwrap(std::vector<R>& vec) {
  auto w = CVecRef<decay_t<R>>{};
  for (const auto& v : vec)
    w.emplace_back(std::cref(v));
  return w;
}

//! Takes a begin and end iterators and returns a vector of references to each element
template <class R, class ForwardIt>
auto wrap(ForwardIt begin, ForwardIt end) {
  auto w = VecRef<R>{};
  for (auto it = begin; it != end; ++it)
    w.emplace_back(std::ref(*it));
  return w;
}

//! Takes a begin and end iterators and returns a vector of references to each element
template <class R, class ForwardIt>
auto cwrap(ForwardIt begin, ForwardIt end) {
  auto w = CVecRef<R>{};
  for (auto it = begin; it != end; ++it)
    w.emplace_back(std::cref(*it));
  return w;
}

/*!
 * @brief Wraps references for each parameter in an iterable container that implements begin()/end()
 * @tparam IterableContainer should have begin()/end() either as free or member functions
 * @tparam R element type
 * @param parameters parameters to wrap
 * @return vector of constant references to each element in parameters
 */
template <class R, class IterableContainer>
auto wrap(const IterableContainer& parameters) {
  using std::begin;
  using std::end;
  auto w = CVecRef<R>{};
  for (auto it = begin(parameters); it != end(parameters); ++it)
    w.emplace_back(std::cref(*it));
  return w;
}

/*!
 * @brief Wraps references for each parameter in an iterable container that implements begin()/end()
 * @tparam IterableContainer should have begin()/end() either as free or member functions
 * @tparam R element type
 * @param parameters parameters to wrap
 * @return vector of references to each element in parameters
 */
template <class R, class IterableContainer>
auto wrap(IterableContainer& parameters) {
  using std::begin;
  using std::end;
  auto w = VecRef<R>{};
  for (auto it = begin(parameters); it != end(parameters); ++it)
    w.emplace_back(std::ref(*it));
  return w;
}

/*!
 * @brief Wraps references for each parameter in an iterable container that implements begin()/end()
 * @tparam IterableContainer should have begin()/end() either as free or member functions
 * @tparam R element type
 * @param parameters parameters to wrap
 * @return vector of constant references to each element in parameters
 */
template <class R, class IterableContainer>
auto cwrap(IterableContainer& parameters) {
  using std::begin;
  using std::end;
  auto w = CVecRef<R>{};
  for (auto it = begin(parameters); it != end(parameters); ++it)
    w.emplace_back(std::cref(*it));
  return w;
}

} // namespace molpro::linalg::itsolv
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_WRAP_H
