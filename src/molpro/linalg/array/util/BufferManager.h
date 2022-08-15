#ifndef ITERATIVE_SOLVER_BUFFERMANAGER_H
#define ITERATIVE_SOLVER_BUFFERMANAGER_H

#include "molpro/linalg/array/DistrArrayFile.h" // TODO make this work with generic DistrArray
#include <future>
#include <molpro/linalg/array/DistrArray.h>
#include <vector>

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
  BufferManager(const CVecRef<DistrArrayFile>& arrays, size_t buffer_size = 8192, int number_of_buffers = 2);
//  BufferManager(const std::vector<DistrArray>& arrays, size_t buffer_size = 8192, int number_of_buffers = 2)
//      : BufferManager(cwrap(arrays), buffer_size, number_of_buffers) {}
//  BufferManager(const DistrArray& array, size_t buffer_size = 8192, int number_of_buffers = 2)
//      : BufferManager(CVecRef<DistrArray>{std::cref(array)}, buffer_size, number_of_buffers) {}

  /*!
   * The current buffer size
   * @return
   */
  size_t buffer_size() const;
  /*!
   * The distance between buffers
   * @return
   */
  size_t buffer_stride() const;
  /*!
   * Offset of the current buffer
   * @return
   */
  size_t buffer_offset() const;
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
  [[nodiscard]] Span<value_type> next(bool initial = false);

  size_t m_current_segment = 0; // current chunk
  const std::pair<size_t, size_t>
      m_range; // memory offsets for MPI found via distr_array_disk.distribution().range(...)
  std::vector<DistrArray::value_type> m_buffer;
  std::future<void> m_read_future;
};

} // namespace molpro::linalg::array::util

#endif // ITERATIVE_SOLVER_BUFFERMANAGER_H
