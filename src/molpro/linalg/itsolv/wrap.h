#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_WRAP_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_WRAP_H
#include <algorithm>
#include <functional>
#include <iterator>
#include <limits>
#include <map>
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

/*!
 * @brief Given wrapped references in wparams and a range of original parameters [begin, end), returns
 * indices of paramters that are wrapped in this range.
 * @param wparams wrapped references to subset of params
 * @param begin start of range of original paramaters
 * @param end end of range of original paramaters
 */
template <typename R, typename ForwardIt>
std::vector<size_t> find_ref(const VecRef<R>& wparams, ForwardIt begin, ForwardIt end) {
  auto indices = std::vector<size_t>{};
  for (auto it = begin; it != end; ++it) {
    auto it_found = std::find_if(std::begin(wparams), std::end(wparams),
                                 [&it](const auto& w) { return std::addressof(w.get()) == std::addressof(*it); });
    if (it_found != std::end(wparams))
      indices.emplace_back(std::distance(begin, it));
  }
  return indices;
}

/*!
 * @brief Removes indices from a vector
 * @tparam T value type
 * @tparam I index type
 * @param params container
 * @param indices indices to remove
 * @return
 */
template <typename T, typename I>
auto remove_elements(std::vector<T> params, const std::vector<I>& indices) {
  const auto n = params.size();
  for (size_t i = 0, j = 0; i < n; ++i) {
    if (std::find(begin(indices), end(indices), i) != end(indices)) {
      params.erase(begin(params) + j);
    } else {
      ++j;
    }
  }
  return params;
}

} // namespace molpro::linalg::itsolv
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_WRAP_H
