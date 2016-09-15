#ifndef PARAMETERVECTOR_H
#define PARAMETERVECTOR_H

#include <cstddef>
#include <vector>
#ifndef nullptr
#define nullptr NULL
#endif


namespace IterativeSolver {
  typedef double ParameterScalar;

  class ParameterVector
  {
  public:
      /*!
     * \brief Construct an object from provided data. The data buffer will be assumed to continue
     * to exist whilst this object exists, but may have its contents changed externally.
     * \param buffer The data buffer.
     * \param length The length of the data buffer.
     * \param variance
     * - -1: covariant vector
     * - +1: contravariant vector
     * - 0: self-dual vector
     * The class is expected to check that appropriate combinations of vectors are provided in methods that perform linear algebra functions.
     */
    ParameterVector(ParameterScalar* buffer, size_t length, int variance=0);
      /*!
     * \brief Construct an object without any data. Subsequent assignment of data to the object via the overloaded = operator should result in the data
     * being stored internally, either in memory or on external storage.
     * \param variance
     * - -1: covariant vector
     * - +1: contravariant vector
     * - 0: self-dual vector
     * The class is expected to check that appropriate combinations of vectors are provided in methods that perform linear algebra functions.
     */
    ParameterVector(int variance=0);
    virtual ~ParameterVector();
    /*!
     * \brief Add a constant times another object to this object
     * \param a The factor to multiply.
     * \param other The object to be added to this.
     * \return
     */
    virtual void axpy(ParameterScalar a, const ParameterVector& other);
    /*!
     * \brief Scalar product of two objects.
     * \param other The object to be contracted with this.
     * \return
     */
    virtual ParameterScalar operator*(const ParameterVector& other) const;
    /*!
     * \brief Set the contents of the object to zero.
     */
    virtual void zero();
    /*!
     * \brief Copy from one object to another, allocating storage if needed
     * \param other The source of data.
     * \return
     */
    virtual ParameterVector& operator=(const ParameterVector& other);
    virtual ParameterScalar& operator[](size_t pos) { return this->buff[pos];}
    const virtual ParameterScalar& operator[](size_t pos) const { return this->buff[pos];}
    virtual size_t size() const {return this->buff.size();}
  protected:
    int variance_;
  private:
    size_t length_;
    ParameterScalar* buffer_;
    std::vector<ParameterScalar> buff;
  };

#include <iostream>
 inline std::ostream& operator<<(std::ostream& os, ParameterVector const& pv) {
//	 os << "ParameterVector object:";
//	 for (size_t k=0; k<pv.size(); k++)
//		 os <<" "<< pv[k];
//	 os << std::endl;
	 return os;
 }

  class ParameterVectorSetDodgy : private std::vector<ParameterVector>
    {
    public:
  	  ParameterVectorSetDodgy() : std::vector<ParameterVector>() {};
          using std::vector<ParameterVector>::size;
          using std::vector<ParameterVector>::operator[];
          using std::vector<ParameterVector>::iterator;
          using std::vector<ParameterVector>::const_iterator;
          using std::vector<ParameterVector>::begin;
          using std::vector<ParameterVector>::end;
          using std::vector<ParameterVector>::push_back;
          using std::vector<ParameterVector>::pop_back;
          using std::vector<ParameterVector>::resize;
          using std::vector<ParameterVector>::front;
          using std::vector<ParameterVector>::back;

        std::vector<bool> active;
    };
  class ParameterVectorSet
    {
    public:
  	  ParameterVectorSet() {}
  	  size_t size() const { return this->pvs.size();}
  	  ParameterVector& operator[](size_t pos) { return this->pvs[pos];}
  	  const ParameterVector& operator[](size_t pos) const { return pvs[pos];}
  	  ParameterVector& front() { return pvs.front();}
  	  const ParameterVector& front() const { return pvs.front();}
  	  ParameterVector& back() { return pvs.back();}
  	  const ParameterVector& back() const { return pvs.back();}
  	  void push_back(const ParameterVector& val) { pvs.push_back(val); active.push_back(true);}

  	  std::vector<bool> active;
    private:
        std::vector<ParameterVector> pvs;
    };
//  inline std::ostream& operator<<(std::ostream& os, const ParameterVectorSet& pvs) {
//	  for (size_t k=0; k<pvs.size(); k++)
//		  os << pvs[k];
//	 return os;
//  }

}

#endif // PARAMETERVECTOR_H
