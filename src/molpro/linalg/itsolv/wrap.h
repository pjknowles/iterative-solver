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

//! Takes a begin and end iterators and returns a vector of references to each element
template <class ForwardIt>
auto wrap(ForwardIt begin, ForwardIt end) {
  using T = decay_t<decay_t<decltype(*begin)>>;
  auto w = VecRef<T>{};
  for (auto it = begin; it != end; ++it)
    w.push_back(std::ref(*it));
  return w;
}

//! Takes a begin and end iterators and returns a vector of references to each element
template <class ForwardIt>
auto cwrap(ForwardIt begin, ForwardIt end) {
  using T = decay_t<decay_t<decltype(*begin)>>;
  auto w = CVecRef<T>{};
  for (auto it = begin; it != end; ++it)
    w.push_back(std::cref(*it));
  return w;
}

/*!
 * @brief Wraps references for each parameter in an iterable container that implements begin()/end()
 * @tparam IterableContainer should have begin()/end() either as free or member functions
 * @tparam R element type
 * @param parameters parameters to wrap
 * @return vector of constant references to each element in parameters
 */
template <class IterableContainer>
auto wrap(const IterableContainer& parameters) {
  using T = typename IterableContainer::value_type;
  using std::begin;
  using std::end;
  auto w = CVecRef<decay_t<T>>{};
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
template <class IterableContainer>
auto wrap(IterableContainer& parameters) {
  using T = typename IterableContainer::value_type;
  using std::begin;
  using std::end;
  auto w = VecRef<decay_t<T>>{};
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
template <class IterableContainer>
auto cwrap(IterableContainer& parameters) {
  using T = typename IterableContainer::value_type;
  using std::begin;
  using std::end;
  auto w = CVecRef<decay_t<T>>{};
  for (auto it = begin(parameters); it != end(parameters); ++it)
    w.emplace_back(std::cref(*it));
  return w;
}

/*!
 * @brief Constructs a vector of reference wrappers with provided arguments
 */
template <typename T, typename... S>
auto wrap_arg(T&& arg, S&&... args)
    -> std::enable_if_t<std::conjunction_v<std::is_same<decay_t<T>, decay_t<S>>...>, VecRef<decay_t<T>>> {
  using std::begin;
  using std::end;
  auto w = VecRef<decay_t<T>>{};
  w.emplace_back(std::ref(arg));
  (w.emplace_back(std::ref(args)), ...);
  return w;
}

/*!
 * @brief Constructs a vector of const reference wrappers with provided arguments
 */
template <typename T, typename... S>
auto cwrap_arg(T&& arg, S&&... args)
    -> std::enable_if_t<std::conjunction_v<std::is_same<decay_t<T>, decay_t<S>>...>, CVecRef<decay_t<T>>> {
  using std::begin;
  using std::end;
  auto w = CVecRef<decay_t<T>>{};
  w.emplace_back(std::cref(arg));
  (w.emplace_back(std::cref(args)), ...);
  return w;
}

} // namespace molpro::linalg::itsolv
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_WRAP_H
