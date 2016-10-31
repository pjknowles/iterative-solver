#ifndef CACHEDPARAMETERVECTOR_H
#define CACHEDPARAMETERVECTOR_H

#include "ParameterVector.h"
#include "Storage.h"

namespace IterativeSolver {

  /*!
   * \brief A class that implements ParameterVector with data held on backing store
   */
  class CachedParameterVector : public ParameterVector
  {
  public:
    /*!
   * \brief Construct an object without any data.
   */
    CachedParameterVector(size_t length=0);
    CachedParameterVector(const CachedParameterVector& source);
    ~CachedParameterVector();
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
    CachedParameterVector& operator=(const CachedParameterVector& other);

    // Every child of ParameterVector needs exactly this
    CachedParameterVector* clone() const { return new CachedParameterVector(*this); }

    /*!
     * \brief Specify a cache size for manipulating the data
     * \param length
     */
    void setCacheSize(size_t length);

  private:
    void init();
    /*!
   * \brief The file to hold the data
   */
//    mutable std::fstream m_file;
    mutable Storage* m_file;
    size_t m_size; //!< How much data
    size_t m_cacheSize; //!< cache size for implementing operations
    mutable std::vector<ParameterScalar> m_cache;
    mutable bool m_cacheDirty;
    mutable size_t m_cacheOffset;
    mutable bool m_cacheEmpty;
    void flushCache() const;
    void write(const ParameterScalar * const buffer, size_t length, size_t offset) const;
    void read(ParameterScalar* buffer, size_t length, size_t offset) const;
  public:
    void put(ParameterScalar* const buffer, size_t length, size_t offset);
    void get(ParameterScalar* buffer, size_t length, size_t offset) const;
    ParameterScalar& operator[](size_t pos);
    ParameterScalar& operator[](size_t pos) const;
    size_t size() const {return m_size;}
    std::string str() const;
  };

}
#endif // CACHEDPARAMETERVECTOR_H