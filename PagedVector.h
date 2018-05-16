#ifndef PAGEDVECTOR_H
#define PAGEDVECTOR_H

#include "LinearAlgebra.h"
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string.h>
#include <cstddef>
#include <unistd.h>
#include <fstream>

#ifdef USE_MPI
#include <mpi.h>
#else
using MPI_Comm = int;
constexpr MPI_Comm MPI_COMM_WORLD=0;
constexpr MPI_Comm MPI_COMM_NULL=0;
#endif
extern MPI_Comm MPI_Comm_PagedVector;
MPI_Comm MPI_Comm_PagedVector=MPI_COMM_WORLD;

namespace LinearAlgebra {

 /*!
   * \brief A class that implements LinearAlgebra::vector<scalar> with data optionally held on backing store, and optionally distributed
   * over MPI ranks
   */
 template <class scalar=double  , class Allocator =
          #ifdef MEMORY_MEMORY_H
           memory::allocator<scalar>
          #else
           std::allocator<scalar>
          #endif
>
 class PagedVector : public LinearAlgebra::vector<scalar>
 {
//  static constexpr size_t default_offline_buffer_size=102400; ///< default buffer size if in offline mode
#define default_offline_buffer_size 5000
 public:
  /*!
   * \brief Construct an object without any data.
   */
  PagedVector(size_t length=0, int option=0, MPI_Comm mpi_communicator=MPI_Comm_PagedVector)
   : LinearAlgebra::vector<scalar>(), m_size(length),
     m_mpi_size(mpi_size()), m_mpi_rank(mpi_rank()), m_communicator(mpi_communicator),
     m_replicated(!(LINEARALGEBRA_CLONE_ADVISE_DISTRIBUTED & option)),
     m_segment_offset(m_replicated ? 0 : ((m_size-1) / m_mpi_size + 1) * m_mpi_rank),
     m_segment_length(m_replicated ? m_size : std::min( (m_size-1) / m_mpi_size + 1, m_size-m_segment_offset)),
     m_cache(m_segment_length, (LINEARALGEBRA_CLONE_ADVISE_OFFLINE & option) ? default_offline_buffer_size : m_segment_length)
  {
//   init(option);
//    std::cout <<" option "<<( option)<<std::endl;
//    std::cout <<"LINEARALGEBRA_CLONE_ADVISE_OFFLINE & option "<<(LINEARALGEBRA_CLONE_ADVISE_OFFLINE & option)<<std::endl;
//    std::cout <<"cache preferred length "<<m_cache.preferred_length<<std::endl;
//   std::cout << m_mpi_rank << " in constructor m_segment_length="<<m_segment_length<<", m_segment_offset="<<m_segment_offset<<std::endl;
  }
  PagedVector(const PagedVector& source, int option=0, MPI_Comm mpi_communicator=MPI_Comm_PagedVector)
   : LinearAlgebra::vector<scalar>(), m_size(source.m_size),
     m_mpi_size(mpi_size()), m_mpi_rank(mpi_rank()), m_communicator(mpi_communicator),
     m_replicated(!(LINEARALGEBRA_CLONE_ADVISE_DISTRIBUTED & option)),
     m_segment_offset(m_replicated ? 0 : ((m_size-1) / m_mpi_size + 1) * m_mpi_rank),
     m_segment_length(m_replicated ? m_size : std::min( (m_size-1) / m_mpi_size + 1, m_size-m_segment_offset)),
     m_cache(m_segment_length, (LINEARALGEBRA_CLONE_ADVISE_OFFLINE & option) ? default_offline_buffer_size : m_segment_length)
  {
//   init(option);
//   std::cout << "in copy constructor, before copy, source: "<<source.str()<<std::endl;
//   std::cout << m_mpi_rank << " in copy constructor m_segment_length="<<m_segment_length<<", m_segment_offset="<<m_segment_offset<<std::endl;
//    std::cout <<" option "<<( option)<<std::endl;
//    std::cout <<"LINEARALGEBRA_CLONE_ADVISE_OFFLINE & option "<<(LINEARALGEBRA_CLONE_ADVISE_OFFLINE & option)<<std::endl;
//    std::cout <<"cache preferred length "<<m_cache.preferred_length<<std::endl;
   *this = source;
//   std::cout << "in copy constructor, after copy, source: "<<source.str()<<std::endl;
  }

  virtual ~PagedVector()
  {}

  // Every child of LinearAlgebra::vector<scalar> needs exactly this
  PagedVector* clone(int option=0) const {// std::cout << "in PagedVector clone, option="<<option<<std::endl;
                                           return new PagedVector(*this, option); }

  /*!
     * \brief Specify a cache size for manipulating the data
     * \param length
     */
  void setCacheSize(size_t length) const
  {
   m_cache.move(0,length);
  }
  /*!
     * \brief Whether a full copy of data is replicated on every MPI process
     * \return
     */
  bool replicated() const { return m_replicated;}

 private:
  int mpi_size() {
#ifdef USE_MPI
   int result;
   MPI_Comm_size(m_communicator, &result);
#else
    int result=1;
#endif
    return result;
  }
  int mpi_rank() {
#ifdef USE_MPI
   int result;
   MPI_Comm_rank(m_communicator, &result);
#else
    int result=0;
#endif
    return result;
  }
  size_t seglength(size_t length){
   return length;
  }
  void init(int option)
  {
   m_replicated = true;
#ifdef USE_MPI
//   std::cout << "option="<<option<<std::endl;
    MPI_Comm_size(m_communicator, &m_mpi_size);
    MPI_Comm_rank(m_communicator, &m_mpi_rank);
   if(LINEARALGEBRA_CLONE_ADVISE_DISTRIBUTED & option) m_replicated=false;
#else
    m_mpi_size=1;
    m_mpi_rank=0;
#endif
//    if (m_mpi_rank==0) std::cout << "m_size="<<m_size<< ", m_mpi_size="<<m_mpi_size<<", m_replicated="<<m_replicated<<std::endl;
    if (m_replicated) {
     m_segment_offset=0;
     m_segment_length=m_size;
    } else {
     m_segment_offset = ((m_size-1) / m_mpi_size + 1) * m_mpi_rank;
     m_segment_length = std::min( (m_size-1) / m_mpi_size + 1, m_size-m_segment_offset);
    }
    m_cache.preferred_length = (LINEARALGEBRA_CLONE_ADVISE_OFFLINE & option) ? default_offline_buffer_size : m_segment_length;
    m_cache.move(0);
//   std::cout << "new PagedVector m_size="<<m_size<<", option="<<option<<" cache size="<<m_cache.length<<std::endl;
//    std::cout << "m_mpi_rank="<<m_mpi_rank<< ", m_segment_offset="<<m_segment_offset<<", m_segment_length="<<m_segment_length<<std::endl;
  }

  size_t m_size; //!< How much data
  const MPI_Comm m_communicator; //!< the MPI communicator for distributed data
  int m_mpi_size;
  int m_mpi_rank;
  bool m_replicated; //!< whether a full copy of data is on every MPI process
  size_t m_segment_offset; //!< offset in the overall data object of this process' data
  size_t m_segment_length; //!< length of this process' data

  struct window {
   mutable size_t offset; ///< the offset in mapped data of the first element of the cache window
   mutable size_t length;///< the size of the cache window
   const size_t preferred_length; ///< the default for the size of the cache window
   const size_t datasize; ///< the size of the vector being mapped
   mutable std::vector<scalar> buffer;
   mutable size_t writes; ///< how many scalars written
   mutable size_t reads; ///< how many scalars read
   bool io; ///< whether backing store is needed
   const scalar* begin() const { return buffer.data();}
   const scalar* end() const { return buffer.data()+length;}
   scalar* begin() {return buffer.data();}
   scalar* end() {return buffer.data()+length;}
   mutable bool dirty;
   mutable std::fstream m_file;
   mutable size_t filesize;
   window(size_t datasize, size_t length=default_offline_buffer_size)
    :  datasize(datasize), preferred_length(length), filesize(0), offset(datasize+1), length(0), dirty(false), writes(0), reads(0) {
//    std::cout << "window constructor datasize="<<datasize<<", length="<<length<<std::endl;
    io = this->datasize > preferred_length;
    if (io) {
     char *tmpname=strdup("tmpfileXXXXXX");
     mkstemp(tmpname);
     m_file.open (tmpname, std::ios::out | std::ios::in | std::ios::binary);
     if (!m_file.is_open() ) throw std::runtime_error(std::string("Cannot open cache file ")+tmpname);
     unlink(tmpname);
     free(tmpname);
    }
    else
     buffer.resize(length);
    move(0,length);
//    std::cout << "window constructor ends, filesize="<<filesize<<", length="<<length<<std::endl;
   }

   void move(const size_t offset, size_t length=0) const {
    if (!length) length=preferred_length;
//    std::cout << "move offset="<<offset<<", length="<<length<<", this->length="<<this->length<<", this->offset="<<this->offset<<" this->io="<<this->io<<std::endl;
    if (dirty && this->length && io) {
     m_file.seekp(this->offset*sizeof(scalar));
//          std::cout << "write to "<<this->offset<<":";for (size_t i=0; i<this->length; i++) std::cout <<" "<<buffer[i]; std::cout <<std::endl;
//     std::cout << "write to "<<this->offset<<std::endl;
     m_file.write( (const char*)buffer.data(), this->length*sizeof(scalar));
     writes += this->length;
     if (filesize < this->offset+this->length) filesize = this->offset+this->length;
    }
    this->offset = offset;
    this->length = std::min(length,static_cast<size_t>(datasize-offset));
    if (io) buffer.resize(this->length);
    //    std::cout << "buffer resized to length="<<this->length<<"; offset="<<offset<<", filesize="<<filesize<<std::endl;
    if (std::min(this->length,static_cast<size_t>(filesize-offset)) && io) {
     m_file.seekg(offset*sizeof(scalar));
     m_file.read((char*)buffer.data(),std::min(this->length,static_cast<size_t>(filesize-offset))*sizeof(scalar));
//     std::cout << "read from "<<this->offset<<":";for (size_t i=0; i<this->length; i++) std::cout <<" "<<buffer[i]; std::cout <<std::endl;
//     std::cout << "read from "<<this->offset<<std::endl;
     reads += std::min(this->length,static_cast<size_t>(filesize-offset));
    }
    dirty=false;
//    std::cout << "move done"<<std::endl;
   }

   void ensure(const size_t offset) const {
//    std::cout << "ensure offset="<<offset<<(offset < this->offset || offset >= this->offset+this->length)<<", preferred_length="<<preferred_length<<std::endl;
    if (offset < this->offset || offset >= this->offset+this->length) move(offset,this->preferred_length);
//    std::cout <<"after move, offset="<<this->offset<<", length="<<this->length<<std::endl;
   }

   ~window() {
    move(filesize,0);
    if (io) m_file.close();
   }

   const window& operator++() const {
//    std::cout << "operator++ entry offset="<<offset<<", length="<<length<<std::endl;
    move(offset+length,length);
//    std::cout << "operator++ exit  offset="<<offset<<", length="<<length<<std::endl;
    return *this;
   }
  private:
   window operator++(int) {return this;}
  };
 public:
  window m_cache;

 private:
  void flushCache() {m_cache.ensure(0);}

 public:

  /*!
   * \brief Place the contents
   * \param buffer
   * \param length
   * \param offset
   */
  void put(const scalar* buffer, size_t length, size_t offset)
  {
   // first of all, focus attention only on that part of buffer which appears in [m_segment_offset,m_segment_offset+m_segment_length)
   size_t buffer_offset=0;
   if (offset < m_segment_offset) { buffer_offset = m_segment_offset-offset; offset = m_segment_offset; length -= m_segment_offset-offset;}
   if (offset+length > std::min(m_size,m_segment_offset+m_segment_length)) length = std::min(m_size,m_segment_offset+m_segment_length)-offset;

   // now make relative to this mpi-rank's segment
   offset-=this->m_segment_offset;
   buffer_offset += this->m_segment_offset;

   // first of all, process that part of the data in the initial cache window
   for (size_t k=std::max(offset,m_cache.offset); k<std::min(offset+length,m_cache.offset+m_cache.length); k++) {
    m_cache.buffer[k-m_cache.offset] = buffer[k+buffer_offset];
    m_cache.dirty=true;
    std::cout <<"in initial window m_cache_buffer["<<k-m_cache.offset<<"]=buffer["<<k+buffer_offset<<"]="<<buffer[k+buffer_offset]<<std::endl;
   }

   // next, process the data appearing before the initial cache window
   size_t initial_cache_offset=m_cache.offset; size_t initial_cache_length=m_cache.length;
   for (m_cache.move(offset); m_cache.length && m_cache.offset<initial_cache_offset; ++m_cache) {
    for (size_t k=m_cache.offset; k<m_cache.length&&k<initial_cache_offset; k++)
     m_cache.buffer[k-m_cache.offset] = buffer[k+buffer_offset];
    std::cout <<"processed preceding window"<<std::endl;
    m_cache.dirty=true;
   }

   // finally, process the data appearing after the initial cache window
   for (m_cache.move(initial_cache_offset+initial_cache_length); m_cache.length && m_cache.offset<offset+length; ++m_cache) {
    for (size_t k=m_cache.offset; k<m_cache.length&&k<offset+length; k++)
     m_cache.buffer[k-m_cache.offset] = buffer[k+buffer_offset];
    std::cout <<"processed following window"<<std::endl;
    m_cache.dirty=true;
   }

  }

  void get(scalar* buffer, size_t length, size_t offset) const
  {
   size_t buffer_offset=0;
   if (offset < m_segment_offset) { buffer_offset = m_segment_offset-offset; offset = m_segment_offset; length -= m_segment_offset-offset;}
   if (offset+length > std::min(m_size,m_segment_offset+m_segment_length)) length = std::min(m_size,m_segment_offset+m_segment_length)-offset;
   offset-=this->m_segment_offset;
   size_t off=offset;
   for (m_cache.move(off); off < offset+length && m_cache.length; ++m_cache, off += m_cache.length)
    for (size_t k=0; k<std::min(offset+length-m_cache.offset,m_cache.length); k++)
     buffer[buffer_offset-offset+off+k] = m_cache.buffer[k];
  }


  const scalar& operator[](size_t pos) const
  {
   if (pos >= m_cache.offset+m_segment_offset+m_cache.length || pos < m_cache.offset+m_segment_offset) { // cache not mapping right sector
//    std::cout << "cache miss"<<std::endl;
    if (pos >= m_segment_offset+m_segment_length || pos < m_segment_offset) throw std::logic_error("operator[] finds index out of range");
    m_cache.move(pos-m_segment_offset);
//    std::cout << "cache offset="<<m_cache.offset<<", cache length=" <<m_cache.length<<std::endl;
   }
   return m_cache.buffer[pos-m_cache.offset-m_segment_offset];
  }

  scalar& operator[](size_t pos)
  {
   scalar* result;
   result = &const_cast<scalar&>(static_cast<const PagedVector*>(this)->operator [](pos));
   m_cache.dirty=true;
   return *result;
  }

  size_t size() const {return m_size;}

  std::string str(int verbosity=0, unsigned int columns=UINT_MAX) const {
   std::ostringstream os; os << "PagedVector object:";
   for (size_t k=0; k<m_segment_length; k++)
    os <<" "<< (*this)[m_segment_offset+k];
   os << std::endl;
   return os.str();
  }

  /*!
   * \brief Add a constant times another object to this object
   * \param a The factor to multiply.
   * \param other The object to be added to this.
   * \return
   */
  void axpy(scalar a, const LinearAlgebra::vector<scalar> &other)
  {
   const PagedVector& othe=dynamic_cast <const PagedVector&> (other);
//   std::cout << "PagedVector::axpy this="<<*this<<std::endl;
//   std::cout << "PagedVector::axpy othe="<<*othe<<std::endl;
   if (this->variance() != othe.variance()) throw std::logic_error("mismatching co/contravariance");
   if (this->m_size != m_size) throw std::logic_error("mismatching lengths");
   if (this->m_replicated != othe.m_replicated) throw std::logic_error("mismatching replication status");
   for (m_cache.ensure(0), othe.m_cache.move(0,m_cache.length); m_cache.length; ++m_cache, ++othe.m_cache ) {
     for (size_t i=0; i<m_cache.length; i++)
      m_cache.buffer[i] += a * othe.m_cache.buffer[i];
     m_cache.dirty = true;
   }
  }

  /*!
   * \brief Scalar product of two objects.
   * \param other The object to be contracted with this.
   * \return
   */
  scalar dot(const vector<scalar> &other) const
  {
   const PagedVector& othe=dynamic_cast <const PagedVector&> (other);
   if (this->variance() * othe.variance() < 0) throw std::logic_error("mismatching co/contravariance");
   if (this->m_size != m_size) throw std::logic_error("mismatching lengths");
   if (this->m_replicated != othe.m_replicated) throw std::logic_error("mismatching replication status");
   scalar result=0;
   if (this == &other) {
//    std::cout <<m_mpi_rank<<" dot self="<<std::endl;
    for (m_cache.ensure(0); m_cache.length; ++m_cache) {
//    std::cout <<m_mpi_rank<<" dot self cache length ="<<m_cache.length<<", offset="<<m_cache.offset<<std::endl;
     for (size_t i=0; i<m_cache.length; i++) {
//      std::cout <<m_mpi_rank<< " take i="<<i<<", buffer[i]="<<m_cache.buffer[i]<<std::endl;
      result += m_cache.buffer[i] * m_cache.buffer[i];
     }
    }
   } else {
//    std::cout <<m_mpi_rank<< " dot replication="<<m_replicated<<othe.m_replicated<<std::endl;
//    std::cout <<m_mpi_rank<< " dot m_segment_offset="<<m_segment_offset<<othe.m_segment_offset<<std::endl;
//    std::cout <<m_mpi_rank<< " dot m_segment_length="<<m_segment_length<<othe.m_segment_length<<std::endl;
    for (m_cache.ensure(0), othe.m_cache.move(0,m_cache.length); m_cache.length; ++m_cache, ++othe.m_cache ) {
     for (size_t i=0; i<m_cache.length; i++) {
//      std::cout <<m_mpi_rank<< " take i="<<i<<", buffer[i]="<<m_cache.buffer[i]<<", other="<<othe.m_cache.buffer[i]<<std::endl;
      result += m_cache.buffer[i] * othe.m_cache.buffer[i];
     }
    }
   }
#ifdef USE_MPI
//    std::cout <<m_mpi_rank<<" dot result before reduce="<<result<<std::endl;
   if (!m_replicated) {
    MPI_Allreduce(&result,&result,1,MPI_DOUBLE,MPI_SUM,m_communicator); // FIXME needs attention for non-double
//    std::cout <<m_mpi_rank<<" dot result after reduce="<<result<<std::endl;
   }
#endif
   return result;
  }

  /*!
     * \brief scal Scale the object by a factor.
     * \param a The factor to scale by.
     */
  void scal(scalar a)
  {
   for (m_cache.ensure(0); m_cache.length; ++m_cache) {
     for (size_t i=0; i<m_cache.length; i++)
      m_cache.buffer[i] *= a;
     m_cache.dirty=true;
   }
  }

  /*!
   * \brief Set the contents of the object to zero.
   */
  void zero()
  {
   for (m_cache.ensure(0); m_cache.length; ++m_cache) {
     for (size_t i=0; i<m_cache.length; i++)
      m_cache.buffer[i] = 0;
     m_cache.dirty=true;
   }
  }

  /*!
   * \brief Copy from one object to another, adjusting size if needed.
   * \param other The source of data.
   * \return
   */
  PagedVector& operator=(const PagedVector& other)
  {
   m_size=other.m_size;
   constexpr bool pr=false;//bool pr=m_size<3;
   if (pr) std::cout << "operator= entry, other="<<other<<std::endl;
   this->setVariance(other.variance());
   if (pr) std::cout <<m_mpi_rank<< " operator=() m_segment_offset="<<m_segment_offset<<"="<<other.m_segment_offset<<" for "<<m_segment_length<<"="<<other.m_segment_length<<std::endl;
   if (pr) std::cout <<m_mpi_rank<< " other="<<other<<std::endl;
   size_t cachelength = std::min(m_cache.preferred_length,other.m_cache.preferred_length);
   size_t off=0, otheroff=0;
   //     m_cache.ensure( (m_replicated == other.m_replicated) ? 0: other.m_segment_offset); std::cout << "m_cache.length="<<m_cache.length<<std::endl;
   //     other.m_cache.ensure( (m_replicated == other.m_replicated) ? 0: m_segment_offset); std::cout << "other.m_cache.length="<<other.m_cache.length<<std::endl;
   m_cache.move( (m_replicated == other.m_replicated) ? 0: other.m_segment_offset,m_cache.io ? cachelength : m_size);
   other.m_cache.move((m_replicated == other.m_replicated) ? 0: m_segment_offset,other.m_cache.io ? cachelength : other.m_size);
   while(m_cache.length && other.m_cache.length && off < m_size && otheroff < other.m_size ) {
    if (pr) std::cout <<m_mpi_rank<< " buffer to copy in range "<<m_cache.offset<<"="<<other.m_cache.offset<<" for "<<m_cache.length<<"="<<other.m_cache.length<<std::endl;
//    for (size_t i=0; i<m_cache.length; i++) std::cout <<" "<<other.m_cache.buffer[otheroff+i]; std::cout <<std::endl;
    for (size_t i=0; i<m_cache.length; i++) {
     m_cache.buffer[off+i] = other.m_cache.buffer[otheroff+i];
//     if (pr) std::cout <<m_mpi_rank<<" buffer["<<off+i<<"]="<<m_cache.buffer[off+i]<<std::endl;
    }
    m_cache.dirty = true;
    if (m_cache.io) ++m_cache; else off+=cachelength;
    if (other.m_cache.io) ++other.m_cache; else otheroff+=cachelength;
    if (pr) std::cout << "end of cache loop, off="<<off<<", otheroff="<<otheroff<<std::endl;
   }
   if (pr) std::cout << "operator= after copy loop, other="<<other<<std::endl;
   if (pr) std::cout << "operator= after copy loop, this="<<*this<<std::endl;
#ifdef USE_MPI
   if (m_replicated && ! other.m_replicated) { // replicated <- distributed
    size_t lenseg = ((m_size-1) / m_mpi_size + 1);
    for (int rank=0; rank < m_mpi_size; rank++) {
     size_t off = lenseg * rank;
     for (m_cache.ensure(off); m_cache.length && off < lenseg*(rank+1); off+= m_cache.length, ++m_cache)
      MPI_Bcast(m_cache.buffer.data(),m_cache.length,MPI_DOUBLE,rank,m_communicator); // needs attention for non-double
    }
   }
   //   if (pr) std::cout << "operator= after bcast, other="<<other<<std::endl;
#endif
   return *this;
  }

  bool operator==(const PagedVector& other)
  {
   if (this->variance() != other.variance()) throw std::logic_error("mismatching co/contravariance");
   if (this->m_size != other.m_size) throw std::logic_error("mismatching lengths");
   int diff=0;
    for (m_cache.ensure( (m_replicated == other.m_replicated) ? 0: other.m_segment_offset),
         other.m_cache.move((m_replicated == other.m_replicated) ? 0: m_segment_offset,m_cache.length);
         m_cache.length && other.m_cache.length;
         ++m_cache, ++other.m_cache ) {
     for (size_t i=0; i<m_cache.length; i++) {
      diff = diff || m_cache.buffer[i] == other.m_cache.buffer[i];
     }
    }
#ifdef USE_MPI
   if (!m_replicated && ! other.m_replicated) {
    for (int rank=0; rank < m_mpi_size; rank++) {
      MPI_Allreduce(&diff,&diff,1,MPI_INT,MPI_SUM,m_communicator);
    }
   }
#endif
   return diff;
  }

 };
 template <class scalar>
 inline std::ostream& operator<<(std::ostream& os, PagedVector<scalar> const& obj) { return os << obj.str(); }

 template <class scalar>
 class PagedVectorTest  {
 public:
  PagedVectorTest(size_t n, int option=3) : status(true) {
   PagedVector<scalar> v1(n,option);
//   if (v1.m_mpi_rank==0) std::cout << "PagedVectorTest, n="<<n<<", replicated="<<v1.replicated()<<std::endl;
   std::vector<scalar> vv;
   for (size_t i=0; i<v1.size(); i++)
    vv.push_back(i);
   v1.put(vv.data(),vv.size(),0);
//   std::cout <<v1.m_mpi_rank<< " v1 after assign: "<<v1<<std::endl;
   PagedVector<scalar> v2(v1,option);
   //  std::cout << "v1 after assigning v2: "<<v1<<std::endl;
   //  std::cout << "v2 after assigning v2: "<<v2<<std::endl;
//   std::cout <<v2.m_mpi_rank<< " v2 after assign: "<<v2<<std::endl;
   v1.m_cache.move(0);
   //  std::cout << "v1 "<<v1.str()<<std::endl;
   //  std::cout << "v2 "<<v2.str()<<std::endl;
   double v1v2error = v2.dot(v1)-(n*(n-1)*(2*n-1))/6;
   double v1v1error = 0;//v1.dot(v1)-(n*(n-1)*(2*n-1))/6;
//   std::cout << "v1.v2 error: " <<v1v2error << std::endl;
//   std::cout << "v1.v1 error: " <<v1v1error << std::endl;
   status = status && v1v1error==0 && v1v2error==0;
//   std::cout <<v2.m_mpi_rank<< " v2 after all: "<<v2<<std::endl;
  }
  bool status;

};
}
#endif // PAGEDVECTOR_H
