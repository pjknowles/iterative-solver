#ifndef SIMPLEVECTOR_H
#define SIMPLEVECTOR_H
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif

#include <algorithm>
#include <cassert>
#include <climits>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <map>
#include <numeric>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <unistd.h>
#include <vector>

namespace molpro {
namespace linalg {

/*!
 * \brief A class that implements a vector container that has the following features:
 * - opaque implementation of BLAS including dot(), axpy(), scal()
 * - additional BLAS overloads to combine with a simple sparse vector represented in std::map
 * - efficient import and export of data ranges
 * - read and read-write access to individual elements
 * - additional functions to find the largest elements, and to print the vector
 * \tparam T the type of elements of the vector
 * \tparam Allocator alternative to std::allocator
 */
template <class T = double, class Allocator = std::allocator<T>>
class SimpleArray {
  using container = std::vector<T, Allocator>;
  container m_buffer;

public:
  typedef double scalar_type; // TODO implement this properly from T
  typedef T value_type;
  explicit SimpleArray(size_t length = 0, const T& value = T()) : m_buffer(length, value) {}
  /*!
   * @brief Copy constructor
   * @param source
   * @param option
   */
  SimpleArray<T, Allocator>(const SimpleArray& source, unsigned int option = 0) : m_buffer(source.m_buffer) {}

  /*!
   * \brief Update a range of the object data with the contents of a provided buffer
   * \param buffer
   * \param length
   * \param offset
   */
  void put(const T* buffer, size_t length, size_t offset) { std::copy(buffer, buffer + length, &m_buffer[offset]); }

  /*!
   * \brief Read a range of the object data into a provided buffer
   * @param buffer
   * @param length
   * @param offset
   */
  void get(T* buffer, size_t length, size_t offset) const {
    std::copy(&m_buffer[offset], &m_buffer[offset + length], buffer);
  }

  /*!
   * @brief Return a reference to an element of the data
   * @param pos Offset of the data
   * @return
   */
  const T& operator[](size_t pos) const { return m_buffer[pos]; }

  /*!
   * @brief Return a reference to an element of the data
   * @param pos Offset of the data
   * @return
   */
  T& operator[](size_t pos) { return m_buffer[pos]; }

  /*!
   * @brief Return the number of elements of data
   * @return
   */
  size_t size() const { return m_buffer.size(); }

  /*!
   * \brief Add a constant times another object to this object
   * \param a The factor to multiply.
   * \param other The object to be added to this.
   * \return
   */
  void axpy(scalar_type a, const SimpleArray<T>& other) {
    assert(this->m_buffer.size() == other.m_buffer.size());
    std::transform(other.m_buffer.begin(), other.m_buffer.end(), m_buffer.begin(), m_buffer.begin(),
                   [a](T x, T y) -> T { return y + a * x; });
  }

  /*!
   * \brief Add a constant times a sparse vector to this object
   * \param a The factor to multiply.
   * \param other The object to be added to this.
   * \return
   */
  void axpy(scalar_type a, const std::map<size_t, T>& other) {
    for (const auto& o : other)
      (*this)[o.first] += a * o.second;
  }

  /*!
   * \brief Scalar product of two objects.
   * \param other The object to be contracted with this.
   * \return
   */
  scalar_type dot(const SimpleArray<T>& other) const {
    assert(this->m_buffer.size() == other.m_buffer.size());
    return std::inner_product(m_buffer.begin(), m_buffer.end(), other.m_buffer.begin(), (scalar_type)0);
  }

  /*!
   * \brief Scalar product with a sparse vector
   * \param other The object to be contracted with this.
   * \return
   */
  scalar_type dot(const std::map<size_t, T>& other) const {
    scalar_type result = 0;
    for (const auto& o : other)
      result += o.second * (*this)[o.first];
    return result;
  }

  /*!
   * \brief scal Scale the object by a factor.
   * \param a The factor to scale by. If a is zero, then the current contents of the object are ignored.
   */
  void scal(scalar_type a) {
    if (a != 0)
      std::transform(m_buffer.begin(), m_buffer.end(), m_buffer.begin(), [a](T& x) { return a * x; });
    else
      m_buffer.assign(m_buffer.size(), 0);
  }

  /*!
   * Find the largest values of the object.
   * @param measure A vector of the same size and matching covariancy, with which the largest contributions to the
   * scalar product with *this are selected.
   * @param maximumNumber At most this number of elements are returned.
   * @param threshold Contributions to the scalar product smaller than this are not included.
   * @return index, value pairs. value is the product of the matrix element and the corresponding element of measure.
   *
   */
  std::tuple<std::vector<size_t>, std::vector<T>>
  select(const SimpleArray<T>& measure, const size_t maximumNumber = 1000, const scalar_type threshold = 0) const {
    std::multimap<T, size_t, std::greater<T>> sortlist;
    if (this == &measure) {
      for (size_t i = 0; i < m_buffer.size(); i++) {
        auto test = m_buffer[i] * m_buffer[i];
        if (test > threshold) {
          sortlist.insert(std::make_pair(test, i));
          if (sortlist.size() > maximumNumber)
            sortlist.erase(std::prev(sortlist.end()));
        }
      }
    } else {
      for (size_t i = 0; i < m_buffer.size(); i++) {
        scalar_type test = m_buffer[i] * measure.m_buffer[i];
        if (test < 0)
          test = -test;
        if (test > threshold) {
          sortlist.insert(std::make_pair(test, i));
          if (sortlist.size() > maximumNumber)
            sortlist.erase(std::prev(sortlist.end()));
        }
      }
    }
    while (sortlist.size() > maximumNumber)
      sortlist.erase(std::prev(sortlist.end()));
    std::vector<size_t> indices;
    indices.reserve(sortlist.size());
    std::vector<T> values;
    values.reserve(sortlist.size());
    for (const auto& p : sortlist) {
      indices.push_back(p.second);
      values.push_back(p.first);
    }
    return std::make_tuple(indices, values);
  }

  bool operator==(const SimpleArray& other) { return m_buffer == other.m_buffer; }
};

} // namespace linalg
} // namespace molpro
#endif // SIMPLEVECTOR_H
