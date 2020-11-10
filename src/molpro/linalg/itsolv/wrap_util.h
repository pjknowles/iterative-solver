#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_WRAP_UTIL_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_WRAP_UTIL_H
#include <molpro/linalg/itsolv/wrap.h>

#include <algorithm>
#include <memory>

namespace molpro::linalg::itsolv {
/*!
 * @brief Given wrapped references in wparams and a range of original parameters [begin, end), returns
 * indices of parameters that are wrapped in this range.
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
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_WRAP_UTIL_H
