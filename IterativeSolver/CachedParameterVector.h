#ifndef CACHEDPARAMETERVECTOR_H
#define CACHEDPARAMETERVECTOR_H

#include "LinearAlgebra.h"
#include "Storage.h"
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string.h>

#ifdef USE_MPI
#include <mpi.h>
#endif

namespace LinearAlgebra {

 const static size_t s_cacheEmpty=std::numeric_limits<size_t>::max();

 /*!
   * \brief A class that implements LinearAlgebra::vector<scalar> with data held on backing store
   */
 template <class scalar=double>
 class CachedParameterVector : public LinearAlgebra::vector<scalar>
 {
 public:
  /*!
   * \brief Construct an object without any data.
   */
  CachedParameterVector(size_t length=0) : LinearAlgebra::vector<scalar>(), m_size(length)
  {
   init();
  }
  CachedParameterVector(const CachedParameterVector& source) : LinearAlgebra::vector<scalar>()
  {
   init();
   *this = source;
  }

  virtual ~CachedParameterVector()
  {if (m_file != nullptr) delete m_file;}
  /*!
   * \brief Add a constant times another object to this object
   * \param a The factor to multiply.
   * \param other The object to be added to this.
   * \return
   */
  void axpy(scalar a, const LinearAlgebra::vector<scalar> *other)
  {
   const CachedParameterVector* othe=dynamic_cast <const CachedParameterVector*> (other);
   if (this->variance() != othe->variance()) throw std::logic_error("mismatching co/contravariance");
   if (false) {
    flushCache();
    std::vector<scalar> buffer(m_cacheSize);
    std::vector<scalar> buffero(m_cacheSize);
    for (size_t block=0; block<m_size; block+=buffer.size()) {
     size_t bs=std::min(buffer.size(),m_size-block);
     read(&buffer[0], bs, block);
     othe->read(&buffero[0], bs, block);
     if (true) {
      for (size_t k=0; k<bs; k++) buffer[k] += a*buffero[k];
     }
     else
      for (size_t k=0; k<bs; k++) (*this)[k+block] += a*buffero[k];
    }
   } else {
    if (m_file == nullptr && othe->m_file == nullptr)
     for (size_t k=0; k<m_size; k++) this->m_cache[k] += a*othe->m_cache[k];
    else
     for (size_t k=0; k<m_size; k++) (*this)[k] += a*(*othe)[k];
   }
  }

  /*!
   * \brief Scalar product of two objects.
   * \param other The object to be contracted with this.
   * \return
   */
  scalar dot(const vector<scalar> *other) const
  {
   const CachedParameterVector* othe=dynamic_cast <const CachedParameterVector*> (other);
   if (this->variance() != othe->variance()) throw std::logic_error("mismatching co/contravariance");
   scalar result=0;
   if (false) {
    std::vector<scalar> buffer(m_cacheSize);
    std::vector<scalar> buffero(m_cacheSize);
    for (size_t block=0; block<m_size; block+=buffer.size()) {
     size_t bs=std::min(buffer.size(),m_size-block);
     read(&buffer[0], bs, block);
     othe->read(&buffero[0], bs, block);
     for (size_t k=0; k<bs; k++) result += buffer[k] * buffero[k];
    }
   } else {
    if (m_file == nullptr && othe->m_file == nullptr)
     for (size_t k=0; k<m_size; k++)
      result += m_cache[k] * othe->m_cache[k];
    else
     for (size_t k=0; k<m_size; k++)
      result += (*this)[k] * (*othe)[k];
   }
   return result;
  }

  /*!
     * \brief scal Scale the object by a factor.
     * \param a The factor to scale by.
     */
  void scal(scalar a)
  {
   if (m_file == nullptr)
    for (size_t k=0; k<m_size; k++)
     m_cache[k] *= a;
   else
    for (size_t k=0; k<m_size; k++)
     (*this)[k] *= a;
  }

  /*!
   * \brief Set the contents of the object to zero.
   */
  void zero()
  {
   if (false) {
    flushCache();
    std::vector<scalar> buffer(m_cacheSize);
    for (size_t k=0; k<m_cacheSize; k++) buffer[k] += 0;
    for (size_t block=0; block<m_size; block+=buffer.size()) {
     size_t bs=std::min(buffer.size(),m_size-block);
     write(&buffer[0], bs, block);
    }
   } else {
    for (size_t k=0; k<m_size; k++) (*this)[k] = 0;
   }
  }

  /*!
   * \brief Copy from one object to another, adjusting size if needed.
   * \param other The source of data.
   * \return
   */
  CachedParameterVector& operator=(const CachedParameterVector& other)
  {
   m_size=other.m_size;
   if (false) {
    std::vector<scalar> buffer(m_cacheSize);
    for (size_t block=0; block<m_size; block+=buffer.size()) {
     size_t bs=std::min(buffer.size(),m_size-block);
     other.read(&buffer[0], bs, block);
     write(&buffer[0], bs, block);
    }
   } else {
    (*this)[0] = other[0]; // to ensure cache initialisation
    if (m_file == nullptr && other.m_file == nullptr)
     for (size_t k=1; k<m_size; k++)
      m_cache[k] = other.m_cache[k];
    else
     for (size_t k=1; k<m_size; k++)
      (*this)[k] = other[k];
   }
   this->setVariance(other.variance());
   return *this;
  }


  // Every child of LinearAlgebra::vector<scalar> needs exactly this
  CachedParameterVector* clone(int option=0) const { return new CachedParameterVector(*this); }

  /*!
     * \brief Specify a cache size for manipulating the data
     * \param length
     */
  void setCacheSize(size_t length) const
  {
   flushCache(true);
   m_cacheSize = length;
   m_cache.resize(m_cacheSize);
   m_cacheOffset=s_cacheEmpty;
   m_cacheDirty=false;
   if (m_cacheSize!=0 && m_cacheSize < m_size) m_file = new Storage; // FIXME parallel
  }
  /*!
     * \brief Whether a full copy of data is replicated on every MPI process
     * \return
     */
  bool replicated() const;
  void setReplicated(bool replicated) const;

 private:
  void init()
  {
   m_file = nullptr;
   m_cacheMax=m_cacheOffset=s_cacheEmpty;
   setCacheSize(0);
  }

  /*!
   * \brief The file to hold the data
   */
  bool m_replicated; //!< whether a full copy of data is on every MPI process
  mutable Storage* m_file; //!< backing store. If nullptr, this means that a file is not being used and everything is in m_cache
  size_t m_size; //!< How much data
  mutable size_t m_cacheSize; //!< cache size for implementing operations
  mutable std::vector<scalar> m_cache;
  mutable bool m_cacheDirty;
  mutable size_t m_cacheOffset;
  mutable size_t m_cacheMax;
  void flushCache(bool force=false) const
  {
   if (m_cacheOffset==s_cacheEmpty) return;
   if (force || (m_cacheDirty && m_file != nullptr)) {
    //          std::cout << "flush Cache offset="<<m_cacheOffset<<" ,length="<<std::min(m_cacheSize,(size_t)size()-m_cacheOffset)<<std::endl;
    //          std::cout << "flush buffer begins "<<m_cache[0]<<std::endl;
    write(&m_cache[0],std::min(m_cacheSize,m_size-m_cacheOffset),m_cacheOffset);
    m_cacheDirty=false;
   }
  }
  void write(const scalar * const buffer, size_t length, size_t offset) const
  {
   //  std::cout << "write "<<length<<std::endl;
   if (m_file == nullptr) m_file = new Storage();
   m_file->write((const char*) buffer,length*sizeof(scalar),offset*sizeof(scalar));
  }

  void read(scalar* buffer, size_t length, size_t offset) const
  {
   //  std::cout << "read  "<<length<<std::endl;
   m_file->read((char*) buffer,length*sizeof(scalar),offset*sizeof(scalar));
  }
 public:
  void put(scalar* const buffer, size_t length, size_t offset)
  {
   if (std::max(m_size,length+offset) <= m_cacheSize) { // FIXME parallel
    for (size_t k=0; k<length; k++) m_cache[k+offset] = buffer[k];
   } else
   {
    flushCache();
    write(buffer,length,offset);
    m_cacheOffset=s_cacheEmpty;
   }
   if (length+offset > m_size) m_size = length+offset;
  }

  void get(scalar* buffer, size_t length, size_t offset) const
  {
   if (m_file == nullptr) {
    for (size_t k=0; k<length; k++) buffer[k] = m_cache[k+offset];
   } else
   {
    flushCache();
    read(buffer,length,offset);
   }
  }


  const scalar& operator[](size_t pos) const
  {
   if (pos >= m_cacheMax || pos < m_cacheOffset) { // cache not mapping right sector
    if (m_cacheSize==0) setCacheSize(m_size); // if no setCacheSize() has been issued, then the default is all in memory
    flushCache();
    m_cacheOffset=pos;
    m_cacheMax=std::min(m_cacheOffset+m_cacheSize,m_size);
    size_t l1 = m_file == nullptr ? 0 : m_file->size()/sizeof(scalar)-m_cacheOffset;
    size_t l2 = std::min(m_cacheSize, m_size-m_cacheOffset);
    size_t l=std::min(l1,l2);
    if (l>0) read(&m_cache[0], l, m_cacheOffset);
   }
   return m_cache[pos-m_cacheOffset];
  }

  scalar& operator[](size_t pos)
  {
   scalar* result;
   result = &const_cast<scalar&>(static_cast<const CachedParameterVector*>(this)->operator [](pos));
   m_cacheDirty=true;
   return *result;
  }

  size_t size() const {return m_size;}
  std::string str(int verbosity=0, unsigned int columns=UINT_MAX) const {
   std::ostringstream os; os << "CachedParameterVector object:";
   flushCache();
   //    std::cout << "@ in str, m_cacheDirty="<<m_cacheDirty<<std::endl;
   for (size_t k=0; k<size(); k++) {
    //        std::cout << "k="<<k<<", m_cacheDirty="<<m_cacheDirty<<std::endl;
    os <<" "<< (*this)[k];
   }
   os << std::endl;
   return os.str();
  }

 };

}
#endif // CACHEDPARAMETERVECTOR_H
