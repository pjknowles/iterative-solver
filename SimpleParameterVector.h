#ifndef SIMPLEPARAMETERVECTOR_H
#define SIMPLEPARAMETERVECTOR_H
#include "ParameterVector.h"

namespace IterativeSolver {

class SimpleParameterVector : public ParameterVector
{
public:
    /*!
   * \brief Construct an object without any data.
   */
  SimpleParameterVector(size_t length=0);
  SimpleParameterVector(const SimpleParameterVector &source);
  ~SimpleParameterVector();
  /*!
   * \brief Add a constant times another object to this object
   * \param a The factor to multiply.
   * \param other The object to be added to this.
   * \return
   */
  void axpy(ParameterScalar a, const ParameterVector *other);
  /*!
   * \brief Scalar product of two objects.
   * \param other The object to be contracted with this.
   * \return
   */
  ParameterScalar dot(const ParameterVector *other) const;
  /*!
   * \brief Set the contents of the object to zero.
   */
  void zero();
  /*!
   * \brief Copy from one object to another, adjusting size if needed.
   * \param other The source of data.
   * \return
   */
  SimpleParameterVector& operator=(const SimpleParameterVector& other);
    /*!
   * \brief Record the co/contra-variance status of the object
   * \param variance
   * - -1: covariant vector
   * - +1: contravariant vector
   * - 0: self-dual vector
   * The class is expected to check that appropriate combinations of vectors are provided in methods that perform linear algebra functions.
   */
  void setVariance(int variance=0) {m_variance=variance;}
  int variance() {return m_variance;}
    SimpleParameterVector* clone() const;
    std::string str() const { std::ostringstream os; os << "SimpleParameterVector object:";
                                      for (size_t k=0; k<size(); k++)
                                          os <<" "<< (*this)[k];
                                      os << std::endl;
                                       return os.str();}
protected:
private:
  int m_variance;
  /*!
   * \brief For a simple implementation, just use an STL vector. Classes that inherit are free to do things differently.
   */
  std::vector<ParameterScalar> m_buffer;
public:
  ParameterScalar& operator[](size_t pos);
  const ParameterScalar& operator[](size_t pos) const;
  size_t size() const;
};

inline std::ostream& operator<<(std::ostream& os, SimpleParameterVector const& pv) {
   os << pv.str();
   return os;
}

}

#endif // SIMPLEPARAMETERVECTOR_H
