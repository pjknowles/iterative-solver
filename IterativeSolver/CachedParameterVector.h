#ifndef CACHEDPARAMETERVECTOR_H
#define CACHEDPARAMETERVECTOR_H

#include "LinearAlgebra.h"
#include "Storage.h"
#include <algorithm>
#include <limits>
#ifndef nullptr
#define nullptr NULL
#endif

#ifdef USE_MPI
#include <mpi.h>
#endif

typedef double ParameterScalar;
typedef LinearAlgebra::vector<ParameterScalar> ParameterVector;
namespace LinearAlgebra {

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
    virtual ~CachedParameterVector();
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
     * \brief scal Scale the object by a factor.
     * \param a The factor to scale by.
     */
    void scal(ParameterScalar a);
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
    void setCacheSize(size_t length) const;
    /*!
     * \brief Whether a full copy of data is replicated on every MPI process
     * \return
     */
    bool replicated() const;
    void setReplicated(bool replicated) const;

  private:
    void init();
    /*!
   * \brief The file to hold the data
   */
    bool m_replicated; //!< whether a full copy of data is on every MPI process
    mutable Storage* m_file; //!< backing store. If nullptr, this means that a file is not being used and everything is in m_cache
    size_t m_size; //!< How much data
    mutable size_t m_cacheSize; //!< cache size for implementing operations
    mutable std::vector<ParameterScalar> m_cache;
    mutable bool m_cacheDirty;
    mutable size_t m_cacheOffset;
    mutable size_t m_cacheMax;
    void flushCache(bool force=false) const;
    void write(const ParameterScalar * const buffer, size_t length, size_t offset) const;
    void read(ParameterScalar* buffer, size_t length, size_t offset) const;
  public:
    void put(ParameterScalar* const buffer, size_t length, size_t offset);
    void get(ParameterScalar* buffer, size_t length, size_t offset) const;

    const ParameterScalar& operator[](size_t pos) const
    {
      if (pos >= m_cacheMax || pos < m_cacheOffset) { // cache not mapping right sector
          if (m_cacheSize==0) setCacheSize(m_size); // if no setCacheSize() has been issued, then the default is all in memory
          flushCache();
          m_cacheOffset=pos;
          m_cacheMax=std::min(m_cacheOffset+m_cacheSize,m_size);
          size_t l1 = m_file == nullptr ? 0 : m_file->size()/sizeof(ParameterScalar)-m_cacheOffset;
          size_t l2 = std::min(m_cacheSize, m_size-m_cacheOffset);
          size_t l=std::min(l1,l2);
          if (l>0) read(&m_cache[0], l, m_cacheOffset);
        }
      return m_cache[pos-m_cacheOffset];
    }

    ParameterScalar& operator[](size_t pos)
    {
      ParameterScalar* result;
      result = &const_cast<ParameterScalar&>(static_cast<const CachedParameterVector*>(this)->operator [](pos));
      m_cacheDirty=true;
      return *result;
    }

    size_t size() const {return m_size;}
    std::string str() const;
  };

}
#endif // CACHEDPARAMETERVECTOR_H
