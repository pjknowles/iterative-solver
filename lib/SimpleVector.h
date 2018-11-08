#ifndef SIMPLEVECTOR_H
#define SIMPLEVECTOR_H
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <climits>
#include <cstring>
#include <cstddef>
#include <cstdlib>
#include <unistd.h>
#include <cassert>
#include <fstream>
#include <sstream>
#include <ostream>
#include <map>
#include <vector>
#include <numeric>

#ifndef LINEARALGEBRA_CLONE_ADVISE_OFFLINE
#define LINEARALGEBRA_CLONE_ADVISE_OFFLINE 0x01
#endif
#ifndef LINEARALGEBRA_CLONE_ADVISE_DISTRIBUTED
#define LINEARALGEBRA_CLONE_ADVISE_DISTRIBUTED 0x02
#endif

#ifdef HAVE_MPI_H
#include <mpi.h>
#ifdef MOLPRO
#include "ppidd.h"
#define MPI_COMM_COMPUTE MPI_Comm_f2c(PPIDD_Worker_comm())
#else
#define MPI_COMM_COMPUTE MPI_COMM_WORLD
#endif
#else
#define MPI_COMM_COMPUTE 0
using MPI_Comm = int;
#endif

namespace LinearAlgebra {

/*!
  * \brief A class that implements a vector container that has the following features:
  * - opaque implementation of BLAS including dot(), axpy(), scal()
  * - additional BLAS overloads to combine with a simple sparse vector represented in std::map
  * - efficient import and export of data ranges
  * - read and read-write access to individual elements
  * - additional functions to find the largest elements, and to print the vector
  * \tparam element_t the type of elements of the vector
  * \tparam Allocator alternative to std::allocator
  */
template<class element_t=double,
    class Allocator =
    std::__1::allocator<element_t>
>
class SimpleVector {
  typedef double scalar_type; //TODO implement this properly from element_t
  std::vector<element_t, Allocator> m_buffer;
 public:
  typedef element_t element_type;
  explicit SimpleVector(size_t length = 0, const element_t& value = element_t())
      : m_buffer(length, value) {}
  /*!
   * @brief Copy constructor
   * @param source
   * @param option
   */
  SimpleVector<element_t, Allocator>(const SimpleVector& source, unsigned int option = 0) : m_buffer(source.m_buffer) {}

  /*!
   * \brief Update a range of the object data with the contents of a provided buffer
   * \param buffer
   * \param length
   * \param offset
   */
  void put(const element_t* buffer, size_t length, size_t offset) {
    std::copy(buffer, buffer + length, &m_buffer[offset]);
  }

  /*!
   * \brief Read a range of the object data into a provided buffer
   * @param buffer
   * @param length
   * @param offset
   */
  void get(element_t* buffer, size_t length, size_t offset) const {
    std::copy(&m_buffer[offset], &m_buffer[offset + length], buffer);
  }

  /*!
   * @brief Return a reference to an element of the data
   * @param pos Offset of the data
   * @return
   */
  const element_t& operator[](size_t pos) const {
    return m_buffer[pos];
  }

  /*!
   * @brief Return a reference to an element of the data
   * @param pos Offset of the data
   * @return
   */
  element_t& operator[](size_t pos) {
    return m_buffer[pos];
  }

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
  void axpy(scalar_type a, const SimpleVector<element_t>& other) {
    assert(this->m_buffer.size() == other.m_buffer.size());
    std::transform(other.m_buffer.begin(),
                   other.m_buffer.end(),
                   m_buffer.begin(),
                   m_buffer.begin(),
                   [a](element_t x, element_t y) -> element_t { return y + a * x; });
  }

  /*!
    * \brief Add a constant times a sparse vector to this object
    * \param a The factor to multiply.
    * \param other The object to be added to this.
    * \return
    */
  void axpy(scalar_type a, const std::map<size_t, element_t>& other) {
    for (const auto& o: other)
      (*this)[o.first] += a * o.second;
  }

  /*!
   * \brief Scalar product of two objects.
   * \param other The object to be contracted with this.
   * \return
   */
  scalar_type dot(const SimpleVector<element_t>& other) const {
    assert(this->m_buffer.size() == other.m_buffer.size());
    return std::inner_product(m_buffer.begin(), m_buffer.end(), other.m_buffer.begin(), (scalar_type)0);
  }

  /*!
   * \brief Scalar product with a sparse vector
   * \param other The object to be contracted with this.
   * \return
   */
  scalar_type dot(const std::map<size_t, element_t>& other) const {
    scalar_type result = 0;
    for (const auto& o: other)
      result += o.second * (*this)[o.first];
    return result;
  }

  /*!
     * \brief scal Scale the object by a factor.
     * \param a The factor to scale by. If a is zero, then the current contents of the object are ignored.
     */
  void scal(scalar_type a) {
    if (a != 0)
      std::transform(m_buffer.begin(), m_buffer.end(), m_buffer.begin(), [a](element_t& x) { return a * x; });
    else
      m_buffer.assign(m_buffer.size(), 0);
  }

  /*!
    * Find the largest values of the object.
    * @param measure A vector of the same size and matching covariancy, with which the largest contributions to the scalar
    * product with *this are selected.
    * @param maximumNumber At most this number of elements are returned.
    * @param threshold Contributions to the scalar product smaller than this are not included.
    * @return index, value pairs. value is the product of the matrix element and the corresponding element of measure.
    *
    */
  std::tuple<std::vector<size_t>, std::vector<element_t> > select(
      const SimpleVector<element_t>& measure,
      const size_t maximumNumber = 1000,
      const scalar_type threshold = 0
  ) const {
    std::multimap<element_t, size_t, std::greater<element_t> > sortlist;
    if (this->m_replicated != measure.m_replicated) throw std::logic_error("mismatching replication status");
    if (this == &measure) {
      for (size_t i = 0; i < m_buffer.size(); i++) {
        auto test = m_buffer[i] * m_buffer[i];
        if (test > threshold) {
          sortlist.insert(std::make_pair(test, i));
          if (sortlist.size() > maximumNumber) sortlist.erase(std::prev(sortlist.end()));
        }
      }
    } else {
      for (size_t i = 0; i < m_buffer.size(); i++) {
        scalar_type test = m_buffer[i] * measure.m_buffer[i];
        if (test < 0) test = -test;
        if (test > threshold) {
          sortlist.insert(std::make_pair(test, i));
          if (sortlist.size() > maximumNumber) sortlist.erase(std::prev(sortlist.end()));
        }
      }
    }
    while (sortlist.size() > maximumNumber) sortlist.erase(std::prev(sortlist.end()));
    std::vector<size_t> indices;
    indices.reserve(sortlist.size());
    std::vector<element_t> values;
    values.reserve(sortlist.size());
    for (const auto& p : sortlist) {
      indices.push_back(p.second);
      values.push_back(p.first);
    }
    return std::make_tuple(indices, values);

  }

  bool operator==(const SimpleVector& other) {
    return m_buffer == other.m_buffer;
  }


};

}
#endif // SIMPLEVECTOR_H
