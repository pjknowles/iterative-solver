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
#else
using MPI_Comm = int;
constexpr MPI_Comm MPI_COMM_WORLD=0;
constexpr MPI_Comm MPI_COMM_NULL=0;
#endif

namespace LinearAlgebra {

 const static size_t s_cacheEmpty=std::numeric_limits<size_t>::max();

 /*!
   * \brief A class that implements LinearAlgebra::vector<scalar> with data optionally held on backing store, and optionally distributed
   * over MPI ranks
   */
 template <class scalar=double , MPI_Comm mpi_communicator = MPI_COMM_WORLD
           , class Allocator =
          #ifdef MEMORY_MEMORY_H
           memory::allocator<scalar>
          #else
           std::allocator<scalar>
          #endif
>
 class CachedParameterVector : public LinearAlgebra::vector<scalar>
 {
 public:
  /*!
   * \brief Construct an object without any data.
   */
  CachedParameterVector(size_t length=0, int option=0) : LinearAlgebra::vector<scalar>(), m_size(length), m_communicator(mpi_communicator)
  {
   init(option);
  }
  CachedParameterVector(const CachedParameterVector& source, int option=0) : LinearAlgebra::vector<scalar>(), m_size(source.m_size), m_communicator(mpi_communicator)
  {
   init(option);
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
    std::vector<scalar,Allocator> buffer(m_cacheSize);
    std::vector<scalar,Allocator> buffero(m_cacheSize);
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
    std::vector<scalar,Allocator> buffer(m_cacheSize);
    std::vector<scalar,Allocator> buffero(m_cacheSize);
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
    std::vector<scalar,Allocator> buffer(m_cacheSize);
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
    std::vector<scalar,Allocator> buffer(m_cacheSize);
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
  CachedParameterVector* clone(int option=0) const { std::cout << "in CachedParameterVector clone, option="<<option<<std::endl;return new CachedParameterVector(*this, option); }

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
  void distribute(MPI_Comm communicator = mpi_communicator) {
   m_communicator = communicator;
#ifdef USE_MPI
#endif
  }

 private:
//  static constexpr size_t default_offline_buffer_size=102400; ///< default buffer size if in offline mode
#define default_offline_buffer_size 1
  void init(int option)
  {
   m_file = nullptr;
   m_cacheMax=m_cacheOffset=s_cacheEmpty;
   if (LINEARALGEBRA_CLONE_ADVISE_OFFLINE & option)
    m_cacheSize=(size_t)default_offline_buffer_size;
#ifdef USE_MPI
   m_replicated = ! LINEARALGEBRA_CLONE_ADVISE_DISTRIBUTED & option;
    MPI_Comm_Size(mpi_communicator, m_mpi_size);
    MPI_Comm_Rank(mpi_communicator, m_mpi_rank);
#else
    m_replicated=true;
    m_mpi_size=1;
    m_mpi_rank=1;
#endif
    setCacheSize(0);
   xout << "new CachedParameterVector m_size="<<m_size<<", option="<<option<<" cache size="<<m_cacheSize<<std::endl;
  }

  bool m_replicated; //!< whether a full copy of data is on every MPI process
  const MPI_Comm m_communicator; //!< the MPI communicator for distributed data
  int m_mpi_size;
  int m_mpi_rank;
  mutable Storage* m_file; //!< backing store. If nullptr, this means that a file is not being used and everything is in m_cache
  size_t m_size; //!< How much data
  mutable size_t m_cacheSize; //!< cache size for implementing operations
  mutable std::vector<scalar,Allocator> m_cache;
  mutable bool m_cacheDirty;
  mutable size_t m_cacheOffset;
  mutable size_t m_cacheMax;

  struct cache {
   const Storage& file;
   const off_t offset;
   const size_t length;
   std::vector<scalar> buffer;
   const scalar* begin() const { return buffer.data();}
   const scalar* end() const { return buffer.data()+length;}
   scalar* begin() {return buffer.data();}
   scalar* end() {return buffer.data()+length;}
   bool dirty;
   mutable std::fstream m_file;
   cache(const Storage& file, const off_t offset, const size_t length)
    : file(file), offset(offset), length(length) {
    dirty=false;
    if (std::min(length,static_cast<size_t>(file.size()-offset)))
     file.read(buffer.data(),std::min(length,static_cast<size_t>(file.size()-offset)),offset);
   }
   ~cache() {
    if (dirty && length)
     file.write(buffer.data(),length,offset);
   }
   cache& operator++() {
    this->~cache();
    this->cache(file,offset+length,length);
    return *this;}
   cache operator++(int) {
    return cache (file,offset+length,length);
  }
  };


  struct cache& cache_begin() { return cache(0,m_cacheSize); }
  struct cache& cache_end() { return cache(m_cacheSize,m_cacheSize); }

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
    xout << "cache miss, m_cacheMax="<<m_cacheMax<<", m_cacheOffset="<<m_cacheOffset <<std::endl;
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
