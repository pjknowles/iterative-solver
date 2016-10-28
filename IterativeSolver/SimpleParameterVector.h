#ifndef SIMPLEPARAMETERVECTOR_H
#define SIMPLEPARAMETERVECTOR_H
#include "ParameterVector.h"

namespace IterativeSolver {

  /*!
   * \brief An implementation of ParameterVector that uses std::vector to store all data.
   */
  class SimpleParameterVector : public ParameterVector
  {
  public:
    /*!
   * \brief Construct an object without any data.
   */
    SimpleParameterVector(size_t length=0);
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

    // Every child of ParameterVector needs exactly this
    SimpleParameterVector* clone() const { return new SimpleParameterVector(*this); }


  private:
    /*!
   * \brief For a simple implementation, just use an STL vector.
   */
    std::vector<ParameterScalar> m_buffer;
  public:
    void put(ParameterScalar* const buffer, size_t length, size_t offset);
    void get(ParameterScalar* buffer, size_t length, size_t offset) const;
    size_t size() const;
    std::string str() const;
  };

}

#endif // SIMPLEPARAMETERVECTOR_H
