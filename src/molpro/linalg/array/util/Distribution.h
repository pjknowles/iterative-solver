#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_DISTRIBUTION_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_DISTRIBUTION_H
#include <algorithm>
#include <cassert>
#include <utility>

namespace molpro {
namespace linalg {
namespace array {
namespace util {

/*!
 * @brief Specifies distribution of a contiguous array into non-overlapping chunks
 * @tparam Ind
 */
template <typename Ind = size_t>
struct Distribution {
  using index_type = Ind;

  Distribution() = default;
  Distribution(const Distribution<Ind> &) = default;
  Distribution(Distribution<Ind> &&) noexcept = default;
  Distribution<Ind> &operator=(const Distribution<Ind> &) = default;
  Distribution<Ind> &operator=(Distribution<Ind> &&) noexcept = default;
  ~Distribution() = default;

  /*!
   * @brief Construct a contiguous distribution list using vector of start indices
   * @param indices list of starting indices for each chunk with last element storing past the end index
   */
  explicit Distribution(const std::vector<index_type> &indices) : chunk_borders{indices} {
    assert(("indices must be sorted", std::is_sorted(begin(chunk_borders, end(chunk_borders)))));
  };

  /*!
   * @brief Maps fist and last index in the array to a pair of chunks encapsulating the corresponding range
   * @param lo first index of range
   * @param hi last index of range
   * @return pair of first and last chunk indices encapsulating this range, or size() if range is outside of border()
   */
  std::pair<int, int> cover(index_type lo, index_type hi) const { return {cover(lo), cover(hi)}; };

  /*!
   * @returns index of chunk containing ind if it is within border or size() otherwise.
   * @param ind array index
   */
  int cover(index_type ind) const {
    if (ind < border().first || ind > border().second)
      return size();
    auto iter = std::upper_bound(begin(chunk_borders), end(chunk_borders), ind);
    int chunk = iter == begin(chunk_borders) ? size() : distance(begin(chunk_borders), prev(iter));
    return chunk;
  };

  //! Returns start and size for section of array local to process_rank
  std::pair<index_type, size_t> range(int chunk) {
    return {chunk_borders[chunk], chunk_borders[chunk + 1] - chunk_borders[chunk]};
  };

  //! Start and end index of the distribution
  std::pair<index_type, index_type> border() { return {chunk_borders.front(), chunk_borders.back()}; }

  //! Number of chunks in the distribution
  size_t size() { return std::max(0, chunk_borders.size() - 1); };

  /*!
   * @brief list of starting indices for each chunk with past the end index as last element.
   * e.g. {|0,1,2,3|,4,5|,6,7,8|9} -> {0,4,6,9,10}
   */
  const std::vector<index_type> chunk_borders;
};

} // namespace util
} // namespace array
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_DISTRIBUTION_H
