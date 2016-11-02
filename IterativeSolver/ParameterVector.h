#ifndef PARAMETERVECTOR_H
#define PARAMETERVECTOR_H

#include <cstddef>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <stdexcept>

namespace LinearAlgebra {
  typedef double ParameterScalar;

  /*!
   * \brief An abstract base class for holding expansion vectors and residual vectors
   * for use by IterativeSolver classes.
   *
   * Concrete implementations should implement all of the public methods, together with ensuring
   * that the
   * copy operator=() performs a deep copy to make a completely independent clone.
   * They must also define a copy method exactly as follows.
   * @code
   *      DerivedParameterVector* clone() const { return new DerivedParameterVector(*this); }
   * @endcode
   * Deriving implementations are free to make their own data storage arrangements, including
   * storing
   * externally and/or distributed across processes, so long as all of the public methods function correctly.
   *
   * The put(), get() and size() methods provided in this base class
   * offer direct access to that data, for which there may be no meaningful implementation,
   * so they should not normally be referenced. However they are defined so that
   * code that uses objects of this class can be tested in conjunction with a derivative
   * class that does provide real implementations of these methods.
   */
  class ParameterVector
  {
  public:
    ParameterVector() : m_variance(0){}
    ParameterVector(const ParameterVector &source){ *this = source; }
    virtual ~ParameterVector(){}
    /*!
     * \brief Add a constant times another object to this object
     * \param a The factor to multiply.
     * \param other The object to be added to this.
     * \return
     */
    virtual void axpy(ParameterScalar a, const ParameterVector* other)=0;
    /*!
     * \brief Scalar product of two objects.
     * \param other The object to be contracted with this.
     * \return
     */
    virtual ParameterScalar dot(const ParameterVector* other) const=0;
    /*!
     * \brief scal Scale the object by a factor.
     * \param a The factor to scale by.
     */
    virtual void scal(ParameterScalar a)=0;
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
    void setVariance(int variance=0) {m_variance=variance;}
    /*!
     * \brief Report the co/contra-variance status of the object
     * \return
     * - -1: covariant vector
     * - +1: contravariant vector
     * - 0: self-dual vector
     */
    int variance() const {return m_variance;}
    /*!
     * \brief Make a copy of the object
     * \return A pointer to the new object
     */
    virtual ParameterVector* clone() const=0;
    /*!
     * \brief Make a printable representation of the object
     * (optional implementation).
     * \return
     */
    virtual std::string str() const {return "";}
    /*!
     * \brief Set some of the object's data.
     * (optional implementation).
     * \param buffer Buffer to provide the data.
     * \param length How many numbers to provide.
     * \param offset Offset of first number.
     * \return
     */
    virtual void put(ParameterScalar* const buffer, size_t length, size_t offset)
    {
      throw std::logic_error("Unimplemented function");
    }
    /*!
     * \brief Fetch some of the object's data.
     * (optional implementation).
     * \param buffer Buffer to receive the data.
     * \param length How many numbers to receive.
     * \param offset Offset of first number.
     * \return
     */
    virtual void get(ParameterScalar* buffer, size_t length, size_t offset) const
    {
      throw std::logic_error("Unimplemented function");
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
  protected:
  private:
    int m_variance;
  };

  /*!
  * \brief Stream a ParameterVector object using its str() method.
  * \param os
  * \param pv
  * \return
  */
  inline std::ostream& operator<<(std::ostream& os, const ParameterVector* pv) { os << pv->str(); return os; }
  /*!
   * \brief A container for a collection of ParameterVector objects
   */
  class ParameterVectorSet
  {
  public:
    ParameterVectorSet() {}
    ParameterVectorSet(const ParameterVectorSet& source)
    {
      for (size_t k=0; k<source.size(); k++) {
          push_back_clone(source.m_pvs[k]);
          m_active[k] = source.m_active[k];
        }
    }

    ~ParameterVectorSet() { while(m_pvs.size()>0) m_pvs.pop_back();}
    size_t size() const { return this->m_pvs.size();}
    ParameterVector* operator[](size_t pos) { return m_pvs[pos];}
    const ParameterVector* operator[](size_t pos) const { return m_pvs[pos];}
    ParameterVector* front() { return (m_pvs.front());}
    const ParameterVector* front() const { return (m_pvs.front());}
    ParameterVector* back() { return (m_pvs.back());}
    const ParameterVector* back() const { return (m_pvs.back());}
    void push_back(ParameterVector* val)
    {
      m_pvs.push_back(val);
      m_owned.push_back(false);
      m_active.push_back(true);
    }

    void push_back_clone(ParameterVector* val)
    {
      m_pvs.push_back(val->clone());
      m_owned.push_back(true);
      m_active.push_back(true);
    }
    void pop_back() { if (m_pvs.size() <=0) return; if (m_owned.back()) delete m_pvs.back(); m_pvs.pop_back(); m_active.pop_back(); m_owned.pop_back();}
    void resize(size_t length) { m_pvs.resize(length); m_active.resize(length);}
    ParameterVectorSet& operator=(const ParameterVectorSet& source)
    {
      while(m_pvs.size()>0) m_pvs.pop_back();
      for (size_t k=0; k<source.size(); k++) {
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
    void axpy(ParameterScalar a, const ParameterVectorSet& other) {
      for (size_t k=0; k<size(); k++) this->m_pvs[k]->axpy(a,other[k]);
    }
    /*!
     * \brief Set the contents of the object to zero.
     */
    void zero() { for (size_t k=0; k<size(); k++) this->m_pvs[k]->zero(); }

    /*!
       * \brief A place to record whether each member of the set is meant to be considered active
       */
    std::vector<bool> m_active;
  private:
    std::vector<bool> m_owned; ///< whether the ParameterVector objects pointed to are owned by this container (and therefore destroyed in the destructor)
    std::vector<ParameterVector*> m_pvs; // contain pointers to support polymorphism
  };

  inline std::ostream& operator<<(std::ostream& os, const ParameterVectorSet& pvs) {
    for (size_t k=0; k<pvs.size(); k++) {
        if (pvs.m_active[k])
          os << pvs[k];
      }
    return os;
  }

}

#endif // PARAMETERVECTOR_H
