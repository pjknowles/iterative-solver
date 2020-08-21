#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_SELECT_MAX_DOT_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_SELECT_MAX_DOT_H
#include <queue>

namespace molpro {
namespace linalg {
namespace array {
namespace util {

/*!
 * @brief Select n indices with largest by absolute value contributions to the dot product
 *
 * @tparam X type of left array, must be iterable complete array
 * @tparam Y type of right array, must be iterable complete array
 * @tparam value_type type for product of X and Y elements
 * @tparam value_type_abs type for absolute value of value_type
 * @param n number of indices to select
 * @param x left array
 * @param y right array
 * @return map of indices and corresponding x,y product
 */
template <class X, class Y, typename value_type, typename value_type_abs>
auto select_max_dot(size_t n, const X& x, const Y& y) {
  using std::abs;
  using std::begin;
  using std::end;
  using std::greater;
  using select_pair = std::pair<value_type_abs, size_t>; // value and index
  auto selection = std::priority_queue<select_pair, std::vector<select_pair>, greater<select_pair>>();
  auto ix = begin(x);
  auto iy = begin(y);
  for (size_t i = 0; i < n; ++i, ++ix, ++iy) {
    selection.emplace(abs((*ix) * (*iy)), i);
  }
  for (size_t i = n; i < std::min(x.size(), y.size()); ++i, ++ix, ++iy) {
    selection.emplace(abs((*ix) * (*iy)), i);
    selection.pop();
  }
  auto selection_map = std::map<size_t, value_type_abs>();
  auto m = selection.size();
  for (size_t i = 0; i < m; ++i) {
    selection_map.emplace(selection.top().second, selection.top().first);
    selection.pop();
  }
  return selection_map;
}

} // namespace util
} // namespace array
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_SELECT_MAX_DOT_H
