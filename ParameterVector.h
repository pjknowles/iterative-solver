#ifndef PARAMETERVECTOR_H
#define PARAMETERVECTOR_H

#include <cstddef>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>

namespace IterativeSolver {
  typedef double ParameterScalar;

  /*!
   * \brief An abstract base class for holding expansion vectors and residual vectors
   * for use by IterativeSolver classes.
   *
   * Deriving implementations should implement all of the public methods, together with ensuring
   * that the
   * copy operator=() performs a deep copy to make a completely independent clone.
   * Deriving implementations are free to make their own data storage arrangements, including
   * storing
   * externally and/or distributed across processes, so long as all of the public methods function correctly.
   * The str(), operator[]() and size() methods provided in this base class
   * offer direct access to that data, for which there may be no meaningful implementation,
   * so they should not normally be referenced. However they are defined so that
   * code that uses objects of this class can be tested in conjunction with a derivative
   * class that does provide real implementations of these methods.
   */
  class ParameterVector
  {
  public:
    ParameterVector(){}
    ParameterVector(const ParameterVector &source){}
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
    virtual void setVariance(int variance=0)=0;
    /*!
     * \brief Report the co/contra-variance status of the object
     * \return
     * - -1: covariant vector
     * - +1: contravariant vector
     * - 0: self-dual vector
     */
    virtual int variance()=0;
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
    virtual std::string str() const =0;
    /*!
     * \brief Access an element of the object's data
     * (optional implementation).
     * \param pos
     * \return
     */
    virtual ParameterScalar& operator[](size_t pos)=0;
    const virtual ParameterScalar& operator[](size_t pos) const=0;
    /*!
     * \brief Report the size of the object's data
     * (optional implementation).
     * \return
     */
    virtual size_t size() const=0;
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
          push_back_clone(source.pvs[k]);
          active[k] = source.active[k];
        }
    }

    ~ParameterVectorSet() { while(pvs.size()>0) pvs.pop_back();}
    size_t size() const { return this->pvs.size();}
    ParameterVector* operator[](size_t pos) { return pvs[pos];}
    const ParameterVector* operator[](size_t pos) const { return pvs[pos];}
    ParameterVector* front() { return (pvs.front());}
    const ParameterVector* front() const { return (pvs.front());}
    ParameterVector* back() { return (pvs.back());}
    const ParameterVector* back() const { return (pvs.back());}
    void push_back(ParameterVector* val)
    {
      pvs.push_back(val);
      owned.push_back(false);
      active.push_back(true);
    }

    void push_back_clone(ParameterVector* val)
    {
      pvs.push_back(val->clone());
      owned.push_back(true);
      active.push_back(true);
    }
    void pop_back() { if (pvs.size() <=0) return; if (owned.back()) delete pvs.back(); pvs.pop_back(); active.pop_back(); owned.pop_back();}
    void resize(size_t length) { pvs.resize(length); active.resize(length);}
    ParameterVectorSet& operator=(const ParameterVectorSet& source)
    {
      while(pvs.size()>0) pvs.pop_back();
      for (size_t k=0; k<source.size(); k++) {
          push_back_clone(source.pvs[k]);
          active[k] = source.active[k];
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
      for (size_t k=0; k<size(); k++) this->pvs[k]->axpy(a,other[k]);
    }
    /*!
     * \brief Set the contents of the object to zero.
     */
    void zero() { for (size_t k=0; k<size(); k++) this->pvs[k]->zero(); }

    /*!
       * \brief A place to record whether each member of the set is meant to be considered active
       */
    std::vector<bool> active;
  private:
    std::vector<bool> owned; ///< whether the ParameterVector objects pointed to are owned by this container (and therefore destroyed in the destructor)
    std::vector<ParameterVector*> pvs; // contain pointers to support polymorphism
  };

  inline std::ostream& operator<<(std::ostream& os, const ParameterVectorSet& pvs) {
    for (size_t k=0; k<pvs.size(); k++) {
        if (pvs.active[k])
          os << pvs[k];
      }
    return os;
  }

}

#endif // PARAMETERVECTOR_H
