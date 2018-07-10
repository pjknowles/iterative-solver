#ifndef LINEARALGEBRA_H
#define LINEARALGEBRA_H

#include <cstddef>
#include <climits>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <memory>
#include <tuple>
#include <map>

namespace LinearAlgebra {

/*!
 * \brief An abstract base class for holding vectors with associated basic linear algebra operations.
 *
 * Concrete implementations should implement all of the virtual public methods, together with ensuring
 * that the
 * copy operator=() and the copy constructor perform a deep copy to make a completely independent clone.
 * Deriving implementations are free to make their own data storage arrangements, including
 * storing
 * externally and/or distributed across processes, so long as all of the public methods function correctly.
 *
 * Implementations must also define a copy method exactly as follows.
 * @code
 *      DerivedVector* clone(int option=0) const { return new DerivedVector(*this, option);  }
 * @endcode
 *
 * The operator[](), put(), get() and size() methods provided in this base class
 * offer an interface for direct access to the object's data,
 * for which there might be no meaningful implementation,
 * so they should not normally be referenced. However they are defined so that
 * code that uses objects of this class can be tested in conjunction with a derivative
 * class that does provide real implementations of these methods. Similarly, a dummy str()
 * virtual method is provided, which can be accessed through the << operator with the intention
 * of providing a printable representation of the object.
 *
 *  \tparam scalar Type of the elements of the vector.
 */
template<class scalar>
class vector {
 public:
  vector<scalar>() : m_variance(0) {}
  /*!
   * \brief copy constructor
   * \param source
   * \param option  might contain advisory information that can be consulted by the implementation:
   * - \c LINEARALGEBRA_CLONE_ADVISE_OFFLINE Slow offline storage can be used, and the client would prefer to save memory by doing so
   * - \c LINEARALGEBRA_CLONE_ADVISE_DISTRIBUTED MPI-distributed storage can be used, and the client would prefer to save memory by doing so
   */
  vector<scalar>(const vector<scalar> &source, int option = 0) { *this = source; }
  virtual ~vector<scalar>() {}
  /*!
   * \brief Add a constant times another object to this object
   * \param a The factor to multiply.
   * \param other The object to be added to this.
   * \return
   */
  virtual void axpy(scalar a, const vector<scalar> &other)=0;

  /*!
   * \brief Add a constant times a sparse vector to this object
   * \param a The factor to multiply.
   * \param other The object to be added to this.
   * \return
   */
  virtual void axpy(scalar a, const std::map<size_t, scalar> &other)=0;

  /*!
   * \brief Scalar product of two objects.
   * \param other The object to be contracted with this.
   * \return
   */
  virtual scalar dot(const vector<scalar> &other) const =0;
  /*!
   * \brief Scalar product with a sparse vector
   * \param other The object to be contracted with this.
   * \return
   */
  virtual scalar dot(const std::map<size_t, scalar> &other) const =0;

  /*!
   * \brief scal Scale the object by a factor.
   * \param a The factor to scale by.
   */
  virtual void scal(scalar a)=0;
  vector<scalar> *operator*=(scalar a) { return (*this) *= a; }

  /*!
   * Find the largest values of the object.
   * @param measure A vector of the same size and matching covariancy, with which the largest contributions to the scalar
   * product with *this are selected.
   * @param maximumNumber At most this number of elements are returned.
   * @param threshold Contributions to the scalar product smaller than this are not included.
   * @return A std::map giving index, value pairs. value is the product of the matrix element and the corresponding element of measure.
   *
   */
  virtual std::tuple<std::vector<size_t>, std::vector<scalar> > select(const vector<scalar> &measure,
                                                                       const size_t maximumNumber = 1000,
                                                                       const scalar threshold = 0) const = 0;

  /*!
   * \brief Set the contents of the object to zero.
   */
  virtual void zero()=0;

  /*!
   * \brief Record the co/contra-variance status of the object
   * \param variance
   * - -1: covariant vector
   * - +1: contravariant vector
   * - 0: self-dual vector
   * The class is expected to check that appropriate combinations of vectors are provided in methods that perform linear algebra functions.
   */
  void setVariance(int variance = 0) { m_variance = variance; }

  /*!
   * \brief Report the co/contra-variance status of the object
   * \return
   * - -1: covariant vector
   * - +1: contravariant vector
   * - 0: self-dual vector
   */
  int variance() const { return m_variance; }

#define LINEARALGEBRA_CLONE_ADVISE_OFFLINE 0x01
#define LINEARALGEBRA_CLONE_ADVISE_DISTRIBUTED 0x02
  /*!
   * \brief Make a copy of the object
   * \param option One or more of the following bit patterns combined:
   * - \c LINEARALGEBRA_CLONE_ADVISE_OFFLINE Suggests that the copy be stored on disk
   * - \c LINEARALGEBRA_CLONE_ADVISE_DISTRIBUTED Suggests that the copy be stored distributed across MPI ranks
   * The implementation doesn't have to follow the suggestions
   * \return A pointer to the new object
   */
  virtual vector<scalar> *clone(int option = 0) const =0;

  /*!
   * \brief Make a printable representation of the object
   * (optional implementation).
   * \param verbosity How much to print.
   * \param columns Page width.
   * \return
   */
  virtual std::string str(int verbosity = 0, unsigned int columns = UINT_MAX) const { return ""; }

  /*!
   * \brief Set some of the object's data.
   * (optional implementation).
   * \param buffer Buffer to provide the data.
   * \param length How many numbers to provide.
   * \param offset Offset of first number.
   * \return
   */
  virtual void put(const scalar *buffer,
                   size_t length,
                   size_t offset) { throw std::logic_error("Unimplemented function"); }

  /*!
   * \brief Fetch some of the object's data.
   * (optional implementation).
   * \param buffer Buffer to receive the data.
   * \param length How many numbers to receive.
   * \param offset Offset of first number.
   * \return
   */
  virtual void get(scalar *buffer,
                   size_t length,
                   size_t offset) const { throw std::logic_error("Unimplemented function"); }

  virtual const scalar &operator[](size_t pos) const { throw std::logic_error("Unimplemented function"); }
  virtual scalar &operator[](size_t pos) {
    return const_cast<scalar &> (static_cast<const vector *> (this)->operator[](pos));
  }

  /*!
   * \brief Report the size of the object's data
   * (optional implementation).
   * \return
   */
  virtual size_t size() const //{return 0;}
  {
    throw std::logic_error("Unimplemented function");
  }

 private:
  int m_variance;
};

/*!
* \brief Stream a vector object using its str() method.
* \param os
* \param pv
* \return
*/
template<class scalar>
inline std::ostream &operator<<(std::ostream &os, const vector<scalar> *pv) {
  os << pv->str();
  return os;
}

/*!
 * \brief A container for a collection of LinearAlgebra::vector objects.
 *
 * This class does not need to be reimplemented for derivatives of LinearAlgebra::vector,
 * but can simply be used as-is; vectors can be added by giving a pointer to push_back()
 * (pointer to original data stored) or push_back_clone (copy of original data taken, and destroyed with the object).
 */
template<class scalar>
class vectorSet {
 public:
  vectorSet<scalar>() {}
  vectorSet<scalar>(const vectorSet<scalar> &source, int option = 0) {
//     std::cout << "vectorSet<scalar> about to push_back_clone option="<<option<<std::endl;
    for (size_t k = 0; k < source.size(); k++) {
      push_back_clone(source.m_pvs[k], option);
      m_active[k] = source.m_active[k];
    }
  }

  typedef std::shared_ptr<vector<scalar> > pv_t;
//    typedef vector<scalar> * pv_t;
  ~vectorSet<scalar>() { while (m_pvs.size() > 0) m_pvs.pop_back(); }
  size_t size() const { return this->m_pvs.size(); }
  pv_t operator[](size_t pos) { return m_pvs[pos]; }
  const pv_t operator[](size_t pos) const { return m_pvs[pos]; }
  pv_t front() { return (m_pvs.front()); }
  const pv_t front() const { return (m_pvs.front()); }
  pv_t back() { return (m_pvs.back()); }
  const pv_t back() const { return (m_pvs.back()); }
  void push_back(const pv_t &val, int option = 0) {
    m_pvs.push_back(val);
    m_owned.push_back(false);
    m_active.push_back(true);
  }

  void push_back_clone(const pv_t &val, int option = 0) {
//     std::cout << "in push_back_clone, about to call clone method with option "<<option<<std::endl;
    m_pvs.push_back(std::shared_ptr<vector<scalar> >(val->clone(option)));
    m_owned.push_back(true);
    m_active.push_back(true);
  }
  void pop_back() {
    if (m_pvs.size() <= 0) return;
    m_pvs.pop_back();
    m_active.pop_back();
    m_owned.pop_back();
  }
  void resize(size_t length) {
    m_pvs.resize(length);
    m_active.resize(length);
  }
  virtual vectorSet<scalar> &operator=(const vectorSet<scalar> &source) {
    while (m_pvs.size() > 0) m_pvs.pop_back();
    for (size_t k = 0; k < source.size(); k++) {
      push_back_clone(source.m_pvs[k]);
      m_active[k] = source.m_active[k];
    }
    return *this;
  }

  /*!
   * \brief Add a constant times another object to this object
   * \param a The factor to multiply.
   * \param other The object to be added to this.
   * \return
   */
  void axpy(scalar a, const vectorSet<scalar> &other) {
    for (size_t k = 0; k < size(); k++) this->m_pvs[k]->axpy(a, *other[k]);
  }

  /*!
   * \brief Set the contents of the object to zero.
   */
  void zero() { for (size_t k = 0; k < size(); k++) this->m_pvs[k]->zero(); }

  /*!
   * \brief A place to record whether each member of the set is meant to be considered active
   */
  std::vector<bool> m_active;
 private:
  std::vector<bool>
      m_owned; ///< whether the vector objects pointed to are owned by this container (and therefore destroyed in the destructor)
  std::vector<pv_t> m_pvs; // contain pointers to support polymorphism
};

template<class scalar>
inline std::ostream &operator<<(std::ostream &os, const vectorSet<scalar> &pvs) {
  for (size_t k = 0; k < pvs.size(); k++) {
    if (pvs.m_active[k])
      os << pvs[k];
  }
  return os;
}

}

#endif // LINEARALGEBRA_H
