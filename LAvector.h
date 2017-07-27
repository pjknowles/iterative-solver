#ifndef LAVECTOR_H
#define LAVECTOR_H
#include <stddef.h>
#include <stdint.h>
#include <sstream>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <assert.h> // assumed that NDEBUG is set properly for optimised compilation
#include <iostream>
#include <initializer_list>
#include "memory.h"
#include "LinearAlgebra.h"


namespace LinearAlgebra {


/*!
 * @brief An implementation of LinearAlgebra::vector using memory::vector
 *
 * \tparam T Type of the elements of the array.
 *
 */
template<typename T=double> class LAvector : LinearAlgebra::vector<T>, memory::vector<T> {

public:
  virtual ~LAvector() {  }

  /*!
   * \brief Add a constant times another object to this object
   * \param a The factor to multiply.
   * \param other The object to be added to this.
   * \return
   */
  virtual void axpy(T a, const LinearAlgebra::vector<T>* other)
  {
    const vector<T>* othe=dynamic_cast <const vector<T>*> (other);
    if (this->variance() != othe->variance()) throw std::logic_error("mismatching co/contravariance");
    for (size_t k=0; k<length; k++) buffer[k] += a*othe->buffer[k];
  }

  /*!
   * \brief Scalar product of two objects.
   * \param other The object to be contracted with this.
   * \return
   */
  virtual T dot(const LinearAlgebra::vector<T>* other) const
  {
    const vector<T>* othe=dynamic_cast <const vector<T>*> (other);
    if (this->variance() * othe->variance() > 0) throw std::logic_error("mismatching co/contravariance");
    T result=0;
    for (size_t k=0; k<length; k++) result += buffer[k] * othe->buffer[k];
    return result;
  }

  vector<T>* clone() const { return new vector<T>(*this); }

  /*!
   * \brief scal Scale the object by a factor.
   * \param a The factor to scale by.
   */
  virtual void scal(T a)
  {
    for (size_t k=0; k<length; k++)  buffer[k] *= a;
  }

  /*!
   * \brief Set the contents of the object to zero.
   */
  virtual void zero()
  {
    for (size_t k=0; k<length; k++)  buffer[k] = 0;
  }

};

}

#endif // LAVECTOR_H
