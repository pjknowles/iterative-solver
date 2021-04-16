#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_DISTRIBUTION_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_DISTRIBUTION_H
#include <algorithm>
#include <cassert>
#include <numeric>
#include <utility>
#include <vector>

namespace molpro::linalg::array::util {

/*!
 * @brief Specifies distribution of a contiguous array into non-overlapping chunks
 * @tparam Ind
 */
template <typename Ind = size_t>
class Distribution {
public:
  using index_type = Ind;

  Distribution() = default;
  Distribution(const Distribution<Ind> &) = default;
  Distribution(Distribution<Ind> &&) noexcept = default;
  Distribution<Ind> &operator=(const Distribution<Ind> &) = default;
  Distribution<Ind> &operator=(Distribution<Ind> &&) noexcept = default;
  ~Distribution() = default;

  friend void swap(Distribution<index_type> &l, Distribution<index_type> &r) {
    using std::swap;
    swap(l.m_chunk_borders, r.m_chunk_borders);
  }

  /*!
   * @brief Construct a contiguous distribution list using vector of start indices
   * @param indices list of starting indices for each chunk with last element storing past the end index
   */
  Distribution(const std::vector<index_type> &indices) : m_chunk_borders{indices} {
    assert("indices must be sorted" && std::is_sorted(begin(m_chunk_borders), end(m_chunk_borders)));
  };

  /*!
   * @brief Maps first and last index in the array to a pair of chunks encapsulating the corresponding range
   * @param lo first index of range
   * @param hi past-the-end index of range
   * @return pair of first and last chunk indices encapsulating this range, or size() if range is outside of border()
   */
  std::pair<int, int> cover(index_type lo, index_type hi) const {
    assert(lo < hi);
    return {cover(lo), cover(hi - 1)};
  };

  /*!
   * @returns index of chunk containing ind if it is within border or size() otherwise.
   * @param ind array index
   */
  int cover(index_type ind) const {
    if (ind < border().first || ind > border().second)
      return size();
    auto iter = std::upper_bound(begin(m_chunk_borders), end(m_chunk_borders), ind);
    int chunk = iter == begin(m_chunk_borders) ? size() : distance(begin(m_chunk_borders), prev(iter));
    return chunk;
  };

  //! Returns [fist, end) indices for section of array assigned to chunk
  std::pair<index_type, index_type> range(int chunk) const {
    return {m_chunk_borders[chunk], m_chunk_borders[chunk + 1]};
  };

  //! @returns [start, end) indices of the distribution
  std::pair<index_type, index_type> border() const { return {m_chunk_borders.front(), m_chunk_borders.back()}; }

  //! Number of chunks in the distribution
  size_t size() const { return m_chunk_borders.size() ? m_chunk_borders.size() - 1 : 0; };

  const std::vector<index_type> &chunk_borders() const { return m_chunk_borders; }

  //! Checks that other distribution is the same
  bool compatible(const Distribution<Ind> &other) const {
    return size() == other.size() &&
           std::inner_product(begin(m_chunk_borders), end(m_chunk_borders), begin(other.m_chunk_borders), true,
                              std::logical_and<>(), std::equal_to<>());
  }

protected:
  /*!
   * @brief list of starting indices for each chunk with past the end index as last element.
   * e.g. {|0,1,2,3|,4,5|,6,7,8|9} -> {0,4,6,9,10}
   */
  std::vector<index_type> m_chunk_borders;
};

/*!
 * @brief Creates a distribution with first (dimension % n_chunks) chunks larger by 1 element.
 * @tparam Ind integer type for indexing the array
 * @param dimension number of elements in array
 * @param n_chunks number of chunks to subdivide the array into
 * @return
 */
template <typename Ind>
Distribution<Ind> make_distribution_spread_remainder(size_t dimension, int n_chunks) {
  assert(n_chunks > 0);
  auto chunk_borders = std::vector<Ind>{0};
  chunk_borders.reserve(n_chunks + 1);
  auto block = dimension / n_chunks;
  auto n_extra = dimension % n_chunks;
  std::fill_n(std::back_inserter(chunk_borders), n_extra, block + 1);
  std::fill_n(std::back_inserter(chunk_borders), n_chunks - n_extra, block);
  std::partial_sum(begin(chunk_borders), end(chunk_borders), begin(chunk_borders));
  return {chunk_borders};
}

} // namespace molpro::linalg::array::util

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_DISTRIBUTION_H
