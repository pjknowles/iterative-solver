#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_SELECT_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_SELECT_H
#include <cmath>
#include <complex>
#include <cstdlib>
#include <queue>

namespace molpro::linalg::array::util {
// virtual std::map<size_t, value_type_abs> select(size_t n, const AL &x, bool max = false, bool abs = false) = 0;
/*!
 * @brief Select n indices with largest (or smallest) actual (or absolute) value
 *
 * @tparam X type of array, must be iterable complete array
 * @tparam value_type type for X elements
 * @tparam value_type_abs type for absolute value of value_type
 * @param n number of indices to select
 * @param x array to examine
 * @param max If true, select largest values, otherwise smallest
 * @param ignore_sign If true, consider std::abs() of elements
 * @return map of indices and corresponding x,y product
 */
template <class X, typename value_type,
          typename std::enable_if<std::is_compound<value_type>::value, value_type>::type* = nullptr>
auto select(size_t n, const X& x, bool max = false, bool ignore_sign = false) {
  throw std::logic_error("select not implemented for complex types");
  return std::map<size_t, value_type>();
}
template <class X, typename value_type,
          typename std::enable_if<!std::is_compound<value_type>::value, value_type>::type* = nullptr>
auto select(size_t n, const X& x, bool max = false, bool ignore_sign = false) {
  using std::abs;
  using std::begin;
  using std::end;
  using std::greater;
  using select_pair = std::pair<value_type, size_t>; // value and index
  auto selection = std::priority_queue<select_pair, std::vector<select_pair>, greater<select_pair>>();
  auto ix = begin(x);
  for (size_t i = 0; i < n; ++i, ++ix) {
    selection.emplace(max ? (ignore_sign ? abs((*ix)) : (*ix)) : (ignore_sign ? -abs((*ix)) : -(*ix)), i);
  }
  for (size_t i = n; i < x.size(); ++i, ++ix) {
    selection.emplace(max ? (ignore_sign ? abs((*ix)) : (*ix)) : (ignore_sign ? -abs((*ix)) : -(*ix)), i);
    selection.pop();
  }
  auto selection_map = std::map<size_t, value_type>();
  auto m = selection.size();
  for (size_t i = 0; i < m; ++i) {
    selection_map.emplace(selection.top().second, selection.top().first);
    selection.pop();
  }
  if (not max)
    for (auto& e : selection_map)
      e.second = -e.second;
  return selection_map;
}

/*!
 * @brief Select n indices with largest (or smallest) actual (or absolute) value
 *
 * @tparam X type of array, must be iterable complete array
 * @tparam value_type type for X elements
 * @param n number of indices to select
 * @param x array to examine
 * @param max If true, select largest values, otherwise smallest
 * @param ignore_sign If true, consider std::abs() of elements
 * @return map of indices and corresponding x,y product
 */
// template <class X, typename value_type>
// auto select_iter_sparse(size_t n, const X& x, bool max = false, bool ignore_sign = false) {
//  return select<X,value_type>(n, x, max, ignore_sign);
//  using std::abs;
//  using std::begin;
//  using std::end;
//  using std::greater;
//  using select_pair = std::pair<value_type, size_t>; // value and index
//  auto selection = std::priority_queue<select_pair, std::vector<select_pair>, greater<select_pair>>();
//  auto ix = begin(x);
//  for (size_t i = 0; i < n; ++i, ++ix) {
//    if (ix->first < x.size())
//      selection.emplace(max ? (ignore_sign ? abs((ix->second)) : (ix->second)) : (ignore_sign ? -abs((ix->second)) :
//      -(ix->second)), ix->first);
//  }
//  for (auto jx = ix; jx != end(x); ++jx) {
//    if (jx->first < x.size()) {
//      selection.emplace(max ? (ignore_sign ? abs((ix->second)) : (ix->second)) : (ignore_sign ? -abs((ix->second)) :
//      -(ix->second)), jx->first); selection.pop();
//    }
//  }
//  auto selection_map = std::map<size_t, value_type>();
//  auto m = selection.size();
//  for (size_t i = 0; i < m; ++i) {
//    selection_map.emplace(selection.top().second, selection.top().first);
//    selection.pop();
//  }
//  return selection_map;
//}

/*!
 * @brief Select n indices with largest (or smallest) actual (or absolute) value
 *
 * @tparam X type of array, must be iterable complete array
 * @tparam value_type type for X elements
 * @tparam value_type_abs type for absolute value of value_type
 * @param n number of indices to select
 * @param x array to examine
 * @param max If true, select largest values, otherwise smallest
 * @param ignore_sign If true, consider std::abs() of elements
 * @return map of indices and corresponding x,y product
 */
template <class X, typename value_type>
auto select_sparse(size_t n, const X& x, bool max = false, bool ignore_sign = false) {
  using std::abs;
  using std::begin;
  using std::end;
  using std::greater;
  using select_pair = std::pair<value_type, size_t>; // value and index
  auto selection = std::priority_queue<select_pair, std::vector<select_pair>, greater<select_pair>>();
  auto ix = begin(x);
  while (selection.size() < n && ix != end(x)) {
    selection.emplace(max ? (ignore_sign ? abs((ix->second)) : (ix->second))
                          : (ignore_sign ? -abs((ix->second)) : -(ix->second)),
                      ix->first);
    ++ix;
  }
  while (ix != end(x)) {
    selection.emplace(max ? (ignore_sign ? abs((ix->second)) : (ix->second))
                          : (ignore_sign ? -abs((ix->second)) : -(ix->second)),
                      ix->first);
    selection.pop();
    ++ix;
  }
  auto selection_map = std::map<size_t, value_type>();
  auto m = selection.size();
  for (size_t i = 0; i < m; ++i) {
    selection_map.emplace(selection.top().second, selection.top().first);
    selection.pop();
  }
  if (not max)
    for (auto& e : selection_map)
      e.second = -e.second;
  return selection_map;
}
} // namespace molpro::linalg::array::util

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_SELECT_H
