#ifndef ITERATIVE_SOLVER_BUFFERMANAGER_H
#define ITERATIVE_SOLVER_BUFFERMANAGER_H

#include "molpro/linalg/array/DistrArrayFile.h" // TODO make this work with generic DistrArray
#include <future>
#include <molpro/linalg/array/DistrArray.h>
#include <vector>
#include <molpro/linalg/array/util/Distribution.h>

namespace molpro::linalg::array::util {

/**
 * @brief BufferManager provides single-buffered or asynchronous double-buffered read access to the data in a collection
 * of DistrArray objects. At construction, an amount of memory is allocated as a buffer, and the buffer is divided into
 * one or more chunks that are independent windows on the data. Sequential access to the data is provided through
 * iterators.
 */
class BufferManager {
  template <class A>
  using VecRef = std::vector<std::reference_wrapper<A>>;

  template <class A>
  using CVecRef = std::vector<std::reference_wrapper<const A>>;
  auto cwrap(const std::vector<DistrArrayFile>& arrays) {
    CVecRef<DistrArrayFile> w;
    for (auto it = arrays.begin(); it != arrays.end(); ++it)
      w.emplace_back(std::cref(*it));
    return w;
  }

public:
  /**
   * @brief Construct a new Buffer Manager object and allocate memory for the chunks.
   *
   * @param arrays the DistrArrayFile objects that this BufferManager will access data from. They must have the same
   * distribution.
   * @param buffer_size number of array elements in each chunk.
   * @param number_of_buffers how many buffers.
   */
  BufferManager(const CVecRef<DistrArrayFile>& arrays, size_t buffer_size = 8192, int number_of_buffers = 2)
      : m_arrays(arrays), m_buffer_size(buffer_size), m_number_of_buffers(number_of_buffers),
        m_range(arrays.empty() ? std::pair<DistrArray::index_type, DistrArray::index_type>{0, 0}
                               : arrays.front().get().distribution().range(molpro::mpi::rank_global())),
        m_buffer(arrays.size() * number_of_buffers * buffer_size) {
    assert(number_of_buffers > 0 && number_of_buffers < 3);
  }

  //  BufferManager(const std::vector<DistrArray>& arrays, size_t buffer_size = 8192, int number_of_buffers = 2)
  //      : BufferManager(cwrap(arrays), buffer_size, number_of_buffers) {}
  //  BufferManager(const DistrArray& array, size_t buffer_size = 8192, int number_of_buffers = 2)
  //      : BufferManager(CVecRef<DistrArray>{std::cref(array)}, buffer_size, number_of_buffers) {}

  /*!
   * The current buffer size
   * @return
   */
  size_t buffer_size() const { return m_current_buffer_size; }
  /*!
   * The distance between buffers
   * @return
   */
  size_t buffer_stride() const { return m_buffer_size; }
  /*!
   * Offset of the current buffer
   * @return
   */
  size_t buffer_offset() const { return m_current_segment * m_buffer_size; }

protected:
  using value_type = DistrArray::value_type;
  const CVecRef<DistrArrayFile>& m_arrays; // reference to the DistrArray objects data is being accessed from
  const size_t m_buffer_size = 8192;
  const int m_number_of_buffers = 2;
  size_t m_current_buffer_size;

public:
  /**
   * @brief Custom iterator for the BufferManager. This iterator is responsible for loading data into the buffers and
   * providing access to that data.
   */
  struct Iterator {
    // iterator properties
    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = Span<DistrArray::value_type>;
    using pointer = Span<DistrArray::value_type>*;
    using reference = const Span<DistrArray::value_type>&;

    /**
     * @brief Construct a new Iterator object
     *
     * @param manager BufferManager object
     * @param begin set this->m_value to the start of the iteration
     * @param end set this->m_value to the end of the iteration
     */
    Iterator(BufferManager& manager, bool begin = false, bool end = false)
        : m_manager(manager), m_value(end ? value_type(nullptr, 0) : manager.next(begin)) {}

    // iterator operators
    reference operator*() const { return m_value; }

    pointer operator->() { return &m_value; }

    Iterator& operator++() {
      m_value = m_manager.next();
      return *this;
    }

    Iterator operator++(int) {
      Iterator tmp = *this;
      ++(*this);
      return tmp;
    }

    friend bool operator==(const Iterator& a, const Iterator& b) {
      return (a.m_value.size() == b.m_value.size() and b.m_value.size() == 0) // special for end of file
             or a.m_value.data() == b.m_value.data();
    };

    friend bool operator!=(const Iterator& a, const Iterator& b) { return not(a == b); };

  private:
    BufferManager& m_manager; // BufferManager object containing buffers
    value_type m_value;       // iterator value
  };

  [[nodiscard]] Iterator begin() { return Iterator(*this, true); }

  [[nodiscard]] Iterator end() { return Iterator(*this, false, true); }

protected:
  /**
   * @brief Move to the next segment of the arrays
   * @param initial whether this is the first iteration
   * @return the iterator value for the next iteration
   */
  [[nodiscard]] Span<value_type> next(bool initial = false) {
    auto prof = molpro::Profiler::single()->push(std::string{"BufferManager::next("} + (initial ? "T" : "F") + ")");
    if (initial)
      m_current_segment = 0;
    else
      ++m_current_segment;
    const auto buffer_id = m_current_segment % m_number_of_buffers;
    const auto lo = m_range.first + m_current_segment * this->m_buffer_size;
    if (lo >= m_range.second)
      return Span<BufferManager::value_type>(nullptr, 0);
    const auto hi = std::min(lo + this->m_buffer_size, m_range.second);
    m_current_buffer_size = hi - lo;
    {
      if (m_number_of_buffers == 1 or lo == m_range.first) {
        for (size_t iarray = 0; iarray < m_arrays.size(); ++iarray) {
          assert(hi > lo);
          assert(buffer_id * m_buffer_size * m_arrays.size() + iarray * m_buffer_size + hi - lo < m_buffer.size());
          m_arrays[iarray].get().get(
              lo, hi, &m_buffer.at(buffer_id * m_buffer_size * m_arrays.size() + iarray * m_buffer_size));
        }
      } else {
        m_read_future.wait();
      }
    }

    const auto next_buffer_id = (m_current_segment + 1) % m_number_of_buffers;
    const auto next_lo = m_range.first + (m_current_segment + 1) * this->m_buffer_size;
    const auto next_hi = std::min(next_lo + this->m_buffer_size, m_range.second);
    if (m_number_of_buffers > 1 and next_lo < m_range.second) {
      m_read_future = std::async(std::launch::async, [this, next_lo, next_hi, next_buffer_id]() {
        for (size_t iarray = 0; iarray < m_arrays.size(); ++iarray)
          m_arrays[iarray].get().get(
              next_lo, next_hi,
              &m_buffer.at(next_buffer_id * m_buffer_size * m_arrays.size() + iarray * m_buffer_size));
      });
    }

    return Span<value_type>(&m_buffer.at(buffer_id * m_buffer_size * m_arrays.size()), m_buffer_size * m_arrays.size());
  }

  size_t m_current_segment = 0; // current chunk
  const std::pair<size_t, size_t>
      m_range; // memory offsets for MPI found via distr_array_disk.distribution().range(...)
  std::vector<DistrArray::value_type> m_buffer;
  std::future<void> m_read_future;
};

} // namespace molpro::linalg::array::util

#endif // ITERATIVE_SOLVER_BUFFERMANAGER_H
