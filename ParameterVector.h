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
     * \brief Calculate the residual for equations to be solved.
     * THIS IS MAYBE NOT GOOD DESIGN SINCE THE CLASS THEN DEPENDS ON WHAT SORT OF EQUATIONS ARE TO BE SOLVED.
     * \param g On exit, will contain the residual added to any previous contents of g
     */
    virtual void residual(ParameterVector & g) const = 0;
    /*!
     * \brief Add a constant times another object to this object
     * \param a The factor to multiply.
     * \param other The object to be added to this.
     * \return
     */
    virtual void axpy(ParameterScalar a, ParameterVector& other);
    /*!
     * \brief Scalar product of two objects.
     * \param other The object to be contracted with this.
     * \return
     */
    virtual ParameterScalar operator*(const ParameterVector& other) const;
    /*!
     * \brief Copy from one object to another, allocating storage if needed
     * \param other The source of data.
     * \return
     */
    virtual ParameterVector& operator=(const ParameterVector& other);
  private:
    int variance_;
    size_t length_;
    ParameterScalar* buffer_;
  };

  class ParameterVectorSet : private std::vector<ParameterVector>
  {

  };

}

#endif // PARAMETERVECTOR_H
