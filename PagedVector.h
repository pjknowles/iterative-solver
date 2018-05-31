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
#include <map>

#ifdef USE_MPI
#include <mpi.h>
#else
using MPI_Comm = int;
constexpr MPI_Comm MPI_COMM_WORLD=0;
#endif

namespace LinearAlgebra {

 /*!
   * \brief A class that implements LinearAlgebra::vector with data optionally held on backing store, and optionally distributed
   * over MPI ranks
   * \tparam scalar the type of elements of the vector
   * \tparam Allocator allocates memory
   */
 template <class scalar=double  , class Allocator =
          #ifdef MEMORY_MEMORY_H
           memory::allocator<scalar>
          #else
           std::allocator<scalar>
          #endif
>
 class PagedVector : public vector<scalar>
 {
//  static constexpr size_t default_offline_buffer_size=102400; ///< default buffer size if in offline mode
#define default_offline_buffer_size 500
 public:
  /*!
   * \brief Construct an object without any data.
   */
  PagedVector(size_t length=0, int option=0, MPI_Comm mpi_communicator=MPI_COMM_WORLD)
   : vector<scalar>(), m_size(length),
     m_communicator(mpi_communicator), m_mpi_size(mpi_size()), m_mpi_rank(mpi_rank()),
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
  PagedVector(const PagedVector& source, int option=0, MPI_Comm mpi_communicator=MPI_COMM_WORLD)
   : vector<scalar>(), m_size(source.m_size),
     m_communicator(mpi_communicator), m_mpi_size(mpi_size()), m_mpi_rank(mpi_rank()),
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
  PagedVector* clone(int option=0) const override {// std::cout << "in PagedVector clone, option="<<option<<std::endl;
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
   const size_t datasize; ///< the size of the vector being mapped
   const size_t preferred_length; ///< the default for the size of the cache window
   mutable std::vector<scalar> buffer;
   bool io; ///< whether backing store is needed
   const scalar* begin() const { return buffer.data();}
   const scalar* end() const { return buffer.data()+length;}
   scalar* begin() {return buffer.data();}
   scalar* end() {return buffer.data()+length;}
   mutable bool dirty;
   mutable std::fstream m_file;
   mutable size_t filesize;
   mutable size_t writes; ///< how many scalars written
   mutable size_t reads; ///< how many scalars read
   window(size_t datasize, size_t length=default_offline_buffer_size)
    :   offset(datasize+1), length(0), datasize(datasize), preferred_length(length), dirty(false), filesize(0), writes(0), reads(0) {
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
//    std::cout << "window move begins"<<std::endl;
//    std::cout << "window move begins, offset="<<offset<<std::endl;
//    std::cout << "window move begins, offset="<<offset<<", length="<<length<<std::endl;
    if (!io) {
     if (offset >= datasize) {
      this->offset=offset;
      this->length=0;
     }
     else {
      this->length=datasize;
      this->offset=0;
     }
//    std::cout << "window move ends !io, offset="<<offset<<", length="<<length<<std::endl;
     return;
    }
    if (!length) length=preferred_length;
//    std::cout << "move offset="<<offset<<", length="<<length<<", this->length="<<this->length<<", this->offset="<<this->offset<<" this->io="<<this->io<<std::endl;
    if (dirty && this->length) {
     m_file.seekp(this->offset*sizeof(scalar));
//          std::cout << "write to "<<this->offset<<":";for (size_t i=0; i<this->length; i++) std::cout <<" "<<buffer[i]; std::cout <<std::endl;
//     std::cout << "write to "<<this->offset<<std::endl;
     m_file.write( (const char*)buffer.data(), this->length*sizeof(scalar));
     writes += this->length;
     if (filesize < this->offset+this->length) filesize = this->offset+this->length;
    }
    this->offset = offset;
    this->length = std::min(length,static_cast<size_t>(datasize-offset));
    buffer.resize(this->length);
    //    std::cout << "buffer resized to length="<<this->length<<"; offset="<<offset<<", filesize="<<filesize<<std::endl;
    if (std::min(this->length,static_cast<size_t>(filesize-offset))) {
     m_file.seekg(offset*sizeof(scalar));
     m_file.read((char*)buffer.data(),std::min(this->length,static_cast<size_t>(filesize-offset))*sizeof(scalar));
//     std::cout << "read from "<<this->offset<<":";for (size_t i=0; i<this->length; i++) std::cout <<" "<<buffer[i]; std::cout <<std::endl;
//     std::cout << "read from "<<this->offset<<std::endl;
     reads += std::min(this->length,static_cast<size_t>(filesize-offset));
    }
    dirty=false;
//    std::cout << "move done"<<std::endl;
//    std::cout << "window move ends, offset="<<offset<<", length="<<length<<std::endl;
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
  void put(const scalar* buffer, size_t length, size_t offset) override
  {
//   std::cout << "PagedVector::put() length="<<length<<", offset="<<offset<<std::endl;
//   for (size_t k=0; k<length; k++) std::cout << " "<<buffer[k]; std::cout << std::endl;
//   std::cout << "cache from "<<m_cache.offset<<" for "<<m_cache.length<<std::endl;
//   for (size_t k=0; k<m_cache.length; k++) std::cout << " "<<m_cache.buffer[k]; std::cout << std::endl;
   // first of all, focus attention only on that part of buffer which appears in [m_segment_offset,m_segment_offset+m_segment_length)
   size_t buffer_offset=0;
   if (offset < m_segment_offset) { buffer_offset = m_segment_offset-offset; offset = m_segment_offset; length -= m_segment_offset-offset;}
   if (offset+length > std::min(m_size,m_segment_offset+m_segment_length)) length = std::min(m_size,m_segment_offset+m_segment_length)-offset;

   // now make addresses relative to this mpi-rank's segment
   offset-=this->m_segment_offset; // the offset in the segment of the first usable element of buffer
   buffer_offset += this->m_segment_offset; // the offset in buffer of its first usable element

//   std::cout << "adjusted offset="<<offset<<", buffer_offset="<<buffer_offset<<std::endl;

   // first of all, process that part of the data in the initial cache window
   for (size_t k=std::max(offset,m_cache.offset); k<std::min(offset+length,m_cache.offset+m_cache.length); k++) { // k is offset in segment
    m_cache.buffer[k-m_cache.offset] = buffer[k+buffer_offset-offset];
    m_cache.dirty=true;
//    std::cout <<"in initial window, k="<<k<<", m_cache_buffer["<<k-m_cache.offset<<"]=buffer["<<k+buffer_offset-offset<<"]="<<buffer[k+buffer_offset-offset]<<std::endl;
   }

   // next, process the data appearing before the initial cache window
   size_t initial_cache_offset=m_cache.offset; size_t initial_cache_length=m_cache.length;
   for (m_cache.move(offset); m_cache.length && m_cache.offset<initial_cache_offset; ++m_cache) {
    for (size_t k=std::max(offset,m_cache.offset); k<m_cache.length&&k<initial_cache_offset&&k<offset+length; k++) { // k is offset in segment
     m_cache.buffer[k-m_cache.offset] = buffer[k+buffer_offset-offset];
    }
//    std::cout <<"processed preceding window"<<std::endl;
    m_cache.dirty=true;
   }

   // finally, process the data appearing after the initial cache window
//std::cout << "initial_cache_offset="<<initial_cache_offset<<std::endl;
   for (m_cache.move(initial_cache_offset+initial_cache_length); m_cache.length && m_cache.offset<offset+length; ++m_cache) {
    for (size_t k=m_cache.offset; k<m_cache.length&&k<offset+length; k++)
     m_cache.buffer[k-m_cache.offset] = buffer[k+buffer_offset-offset];
//    std::cout <<"processed following window"<<std::endl;
    m_cache.dirty=true;
   }

//   std::cout << "PagedVector::put() ends length="<<length<<", offset="<<offset<<std::endl;
  }

  void get(scalar* buffer, size_t length, size_t offset) const override
  {
//   std::cout << "PagedVector::get() length="<<length<<", offset="<<offset<<std::endl;
//   std::cout << "cache from "<<m_cache.offset<<" for "<<m_cache.length<<std::endl;
//   for (size_t k=0; k<m_cache.length; k++) std::cout << " "<<m_cache.buffer[k]; std::cout << std::endl;

   // first of all, focus attention only on that part of buffer which appears in [m_segment_offset,m_segment_offset+m_segment_length)
   size_t buffer_offset=0;
   if (offset < m_segment_offset) { buffer_offset = m_segment_offset-offset; offset = m_segment_offset; length -= m_segment_offset-offset;}
   if (offset+length > std::min(m_size,m_segment_offset+m_segment_length)) length = std::min(m_size,m_segment_offset+m_segment_length)-offset;

   // now make addresses relative to this mpi-rank's segment
   offset-=this->m_segment_offset; // the offset in the segment of the first usable element of buffer
   buffer_offset += this->m_segment_offset; // the offset in buffer of its first usable element

//   std::cout << "adjusted offset="<<offset<<", buffer_offset="<<buffer_offset<<std::endl;

   // first of all, process that part of the data in the initial cache window
   for (size_t k=std::max(offset,m_cache.offset); k<std::min(offset+length,m_cache.offset+m_cache.length); k++) { // k is offset in segment
    buffer[k+buffer_offset-offset]= m_cache.buffer[k-m_cache.offset];
//    std::cout <<"in initial window, k="<<k<<", m_cache_buffer["<<k-m_cache.offset<<"]=buffer["<<k+buffer_offset-offset<<"]="<<buffer[k+buffer_offset-offset]<<std::endl;
   }

   // next, process the data appearing before the initial cache window
   size_t initial_cache_offset=m_cache.offset; size_t initial_cache_length=m_cache.length;
   for (m_cache.move(offset); m_cache.length && m_cache.offset<initial_cache_offset; ++m_cache) {
//    std::cout <<"new cache window offset="<<m_cache.offset<<", length="<<m_cache.length<<std::endl;
    for (size_t k=std::max(offset,m_cache.offset); k<m_cache.length&&k<initial_cache_offset&&k<offset+length; k++) { // k is offset in segment
//     std::cout << "in loop k="<<k<<", buffer_offset="<<buffer_offset<<", offset="<<offset<<" index="<<k+buffer_offset-offset<<std::endl;
     buffer[k+buffer_offset-offset]= m_cache.buffer[k-m_cache.offset];
//    std::cout <<"in preceding window, k="<<k<<", m_cache_buffer["<<k-m_cache.offset<<"]=buffer["<<k+buffer_offset-offset<<"]="<<buffer[k+buffer_offset-offset]<<std::endl;
    }
//    std::cout <<"processed preceding window"<<std::endl;
   }

   // finally, process the data appearing after the initial cache window
   for (m_cache.move(initial_cache_offset+initial_cache_length); m_cache.length && m_cache.offset<offset+length; ++m_cache) {
    for (size_t k=m_cache.offset; k<m_cache.length&&k<offset+length; k++)
     buffer[k+buffer_offset-offset]= m_cache.buffer[k-m_cache.offset];
//    std::cout <<"processed following window"<<std::endl;
   }

//   for (size_t k=0; k<length; k++) std::cout << " "<<buffer[k]; std::cout << std::endl;
//   std::cout << "PagedVector::get() ends, length="<<length<<", offset="<<offset<<std::endl;
  }

  const scalar& operator[](size_t pos) const override
  {
   if (pos >= m_cache.offset+m_segment_offset+m_cache.length || pos < m_cache.offset+m_segment_offset) { // cache not mapping right sector
//    std::cout << "cache miss"<<std::endl;
    if (pos >= m_segment_offset+m_segment_length || pos < m_segment_offset) throw std::logic_error("operator[] finds index out of range");
    m_cache.move(pos-m_segment_offset);
//    std::cout << "cache offset="<<m_cache.offset<<", cache length=" <<m_cache.length<<std::endl;
   }
   return m_cache.buffer[pos-m_cache.offset-m_segment_offset];
  }

  scalar& operator[](size_t pos) override
  {
   scalar* result;
   result = &const_cast<scalar&>(static_cast<const PagedVector*>(this)->operator [](pos));
   m_cache.dirty=true;
   return *result;
  }

  size_t size() const override {return m_size;}

  std::string str(int verbosity=0, unsigned int columns=UINT_MAX) const override {
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
  void axpy(scalar a, const vector<scalar> &other) override
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
  scalar dot(const vector<scalar> &other) const override
  {
//   std::cout << "enter dot"<<std::endl;
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
    double resultLocal = result;
    MPI_Allreduce(&resultLocal,&result,1,MPI_DOUBLE,MPI_SUM,m_communicator); // FIXME needs attention for non-double
//    std::cout <<m_mpi_rank<<" dot result after reduce="<<result<<std::endl;
   }
#endif
//   std::cout << "leave dot"<<std::endl;
   return result;
  }

  /*!
   * \brief Scalar product with a sparse vector
   * \param other The object to be contracted with this.
   * \return
   */
  scalar dot(const std::map<size_t,scalar> &other) const override
  {
   scalar result=0;
   for (const auto& o: other)
    if (o.first >= m_segment_offset && o.first < m_segment_offset+m_segment_length)
     result += o.second * (*this)[o.first];
#ifdef USE_MPI
   if (!m_replicated) {
    double resultLocal = result;
    MPI_Allreduce(&resultLocal,&result,1,MPI_DOUBLE,MPI_SUM,m_communicator); // FIXME needs attention for non-double
   }
#endif
   return result;
  }

  /*!
     * \brief scal Scale the object by a factor.
     * \param a The factor to scale by.
     */
  void scal(scalar a) override
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
  void zero() override
  {
   for (m_cache.ensure(0); m_cache.length; ++m_cache) {
     for (size_t i=0; i<m_cache.length; i++)
      m_cache.buffer[i] = 0;
     m_cache.dirty=true;
   }
  }

  /*!
    * Find the largest values of the object.
    * @param measure A vector of the same size and matching covariancy, with which the largest contributions to the scalar
    * product with *this are selected.
    * @param maximumNumber At most this number of elements are returned.
    * @param threshold Contributions to the scalar product smaller than this are not included.
    * @return index, value pairs. value is the product of the matrix element and the corresponding element of measure.
    *
    */
  std::tuple<std::vector<size_t>,std::vector<scalar> > select (
    const vector <scalar> &measure,
    const size_t maximumNumber = 1000,
    const scalar threshold = 0
  ) const override
  {
   std::multimap<scalar,size_t,std::greater<scalar> > sortlist;
   const PagedVector& measur =dynamic_cast <const PagedVector&> (measure);
   if (this->variance() * measur.variance() < 0) throw std::logic_error("mismatching co/contravariance");
   if (this->m_size != m_size) throw std::logic_error("mismatching lengths");
   if (this->m_replicated != measur.m_replicated) throw std::logic_error("mismatching replication status");
   if (this == &measur) {
    for (m_cache.ensure(0); m_cache.length; ++m_cache) {
     for (size_t i=0; i<m_cache.length; i++) {
      auto test = m_cache.buffer[i] * m_cache.buffer[i];
      if (test > threshold) {
       sortlist.insert(std::make_pair(test,i));
       if (sortlist.size()>maximumNumber) sortlist.erase(std::prev(sortlist.end()));
      }
     }
    }
   } else {
    for (m_cache.ensure(0), measur.m_cache.move(0,m_cache.length); m_cache.length; ++m_cache, ++measur.m_cache ) {
     for (size_t i=0; i<m_cache.length; i++) {
      scalar test = m_cache.buffer[i] * measur.m_cache.buffer[i];
      if (test<0) test=-test;
      if (test > threshold) {
       sortlist.insert(std::make_pair(test,i));
       if (sortlist.size()>maximumNumber) sortlist.erase(std::prev(sortlist.end()));
      }
     }
    }
   }
#ifdef USE_MPI
   if (!m_replicated) {
    throw std::logic_error("incomplete MPI implementation");
   }
#endif
//   std::cout << "leave dot"<<std::endl;
   while (sortlist.size()>maximumNumber) sortlist.erase(std::prev(sortlist.end()));
   std::vector<size_t> indices; indices.reserve(sortlist.size());
   std::vector<scalar> values; values.reserve(sortlist.size());
   for (const auto& p : sortlist) {
    indices.push_back(p.second);
    values.push_back(p.first);
   }
   return std::make_tuple(indices,values);

  };


  /*!
   * \brief Copy from one object to another, adjusting size if needed.
   * \param other The source of data.
   * \return
   */
  PagedVector& operator=(const PagedVector& other)
  {
   m_size=other.m_size;
   constexpr bool pr=false;//bool pr=m_size<3;
//   bool pr=m_mpi_rank>0;
//   if (pr) std::cout << "operator= entry, other="<<other<<std::endl;
   if (pr) std::cout << "operator= entry, other.m_cache.preferred_length="<<other.m_cache.preferred_length<<", this->m_cache.preferred_length="<<this->m_cache.preferred_length<<std::endl;
   if (pr) std::cout << "operator= entry, other.m_replicated="<<other.m_replicated<<", this->m_replicated="<<this->m_replicated<<std::endl;
   this->setVariance(other.variance());
   if (pr) std::cout <<m_mpi_rank<< " operator=() m_segment_offset="<<m_segment_offset<<"="<<other.m_segment_offset<<" for "<<m_segment_length<<"="<<other.m_segment_length<<std::endl;
//   if (pr) std::cout <<m_mpi_rank<< " other="<<other<<std::endl;
   if (!m_cache.io && !other.m_cache.io) { // std::cout << "both in memory"<<std::endl;
    size_t off = (m_replicated && !other.m_replicated) ? other.m_segment_offset : 0;
    size_t otheroff = (!m_replicated && other.m_replicated) ? m_segment_offset : 0;
     for (size_t i=0; i<std::min(m_segment_length,other.m_segment_length); i++) {
      m_cache.buffer[off+i] = other.m_cache.buffer[otheroff+i];
           if (pr) std::cout <<m_mpi_rank<<" other.buffer["<<otheroff+i<<"]="<<other.m_cache.buffer[otheroff+i]<<std::endl;
           if (pr) std::cout <<m_mpi_rank<<" buffer["<<off+i<<"]="<<m_cache.buffer[off+i]<<std::endl;
     }
   } else if (m_cache.io && other.m_cache.io) { // std::cout << "both cached"<<std::endl;
    size_t cachelength = std::min(m_cache.preferred_length,other.m_cache.preferred_length);
//    std::cout << m_mpi_rank<<"cachelength="<<cachelength<<std::endl;
    size_t off=0, otheroff=0;
    //     m_cache.ensure( (m_replicated == other.m_replicated) ? 0: other.m_segment_offset); std::cout << "m_cache.length="<<m_cache.length<<std::endl;
    //     other.m_cache.ensure( (m_replicated == other.m_replicated) ? 0: m_segment_offset); std::cout << "other.m_cache.length="<<other.m_cache.length<<std::endl;
    if (pr) std::cout << "move cache to "<< ((!m_replicated ||  other.m_replicated) ? 0: other.m_segment_offset)<<std::endl;
    if (pr) std::cout << "move other cache to "<<((m_replicated || !other.m_replicated) ? 0: m_segment_offset)<<std::endl;
    m_cache.move( (!m_replicated ||  other.m_replicated) ? 0: other.m_segment_offset, cachelength);
    other.m_cache.move((m_replicated || !other.m_replicated) ? 0: m_segment_offset, cachelength);
    while(m_cache.length && other.m_cache.length && off < m_segment_length && otheroff < other.m_segment_length ) {
     if (pr) std::cout <<m_mpi_rank<< " buffer to copy in range "<<m_cache.offset<<"="<<other.m_cache.offset<<" for "<<m_cache.length<<"="<<other.m_cache.length<<std::endl;
     if (pr) for (size_t i=0; i<m_cache.length; i++) std::cout <<" "<<other.m_cache.buffer[otheroff+i]; std::cout <<std::endl;
     for (size_t i=0; i<m_cache.length; i++) {
      m_cache.buffer[off+i] = other.m_cache.buffer[otheroff+i];
           if (pr) std::cout <<m_mpi_rank<<" buffer["<<off+i<<"]="<<m_cache.buffer[off+i]<<std::endl;
     }
     m_cache.dirty=true;
     if (m_cache.io) ++m_cache; else off+=cachelength;
     if (other.m_cache.io) ++other.m_cache; else otheroff+=cachelength;
     if (pr) std::cout <<m_mpi_rank<< "end of cache loop, off="<<off<<", otheroff="<<otheroff<<std::endl;
    }
   } else if (m_cache.io) { // std::cout<< "source in memory, result cached"<<std::endl;
    size_t cachelength = m_cache.preferred_length;
    size_t off=0;
    size_t otheroff=(m_replicated == other.m_replicated) ? 0: m_segment_offset;
    m_cache.move( (m_replicated == other.m_replicated) ? 0: other.m_segment_offset, cachelength );
    other.m_cache.move( 0, other.m_segment_length);
    while(m_cache.length && off < m_segment_length && otheroff < other.m_segment_length ) {
//      std::cout <<m_mpi_rank<< " buffer to copy in range "<<m_cache.offset<<"="<<other.m_cache.offset<<" for "<<m_cache.length<<"="<<other.m_cache.length<<std::endl;
//      for (size_t i=0; i<m_cache.length; i++) std::cout <<" "<<other.m_cache.buffer[otheroff+i]; std::cout <<std::endl;
//      std::cout <<"off = "<<off<<std::endl;
     for (size_t i=0; i<m_cache.length; i++) {
      m_cache.buffer[off+i] = other.m_cache.buffer[otheroff+i];
      //     if (pr) std::cout <<m_mpi_rank<<" buffer["<<off+i<<"]="<<m_cache.buffer[off+i]<<std::endl;
     }
      m_cache.dirty=true;
     ++m_cache;
     otheroff+=cachelength;
     if (pr) std::cout << "end of cache loop, off="<<off<<", otheroff="<<otheroff<<std::endl;
    }
   } else if (other.m_cache.io) {  // std::cout << "source cached, result in memory"<<std::endl;
    size_t cachelength = other.m_cache.preferred_length;
    size_t otheroff=0;
    size_t off=(m_replicated == other.m_replicated) ? 0: other.m_segment_offset;
    m_cache.move(0, m_segment_length);
    other.m_cache.move((m_replicated == other.m_replicated) ? 0: m_segment_offset, cachelength );
    while(other.m_cache.length && off < m_segment_length && otheroff < other.m_segment_length ) {
     if (pr) std::cout <<m_mpi_rank<< " buffer to copy in range "<<m_cache.offset<<"="<<other.m_cache.offset<<" for "<<m_cache.length<<"="<<other.m_cache.length<<std::endl;
     //    for (size_t i=0; i<m_cache.length; i++) std::cout <<" "<<other.m_cache.buffer[otheroff+i]; std::cout <<std::endl;
     for (size_t i=0; i<other.m_cache.length; i++) {
      m_cache.buffer[off+i] = other.m_cache.buffer[otheroff+i];
      //     if (pr) std::cout <<m_mpi_rank<<" buffer["<<off+i<<"]="<<m_cache.buffer[off+i]<<std::endl;
     }
     off+=cachelength;
     ++other.m_cache;
     if (pr) std::cout << "end of cache loop, off="<<off<<", otheroff="<<otheroff<<std::endl;
    }
   }
   m_cache.dirty = true;
//   if (pr) std::cout << "operator= after copy loop, other="<<other<<std::endl;
//   if (pr) std::cout << "operator= after copy loop, this="<<*this<<std::endl;
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
//std::cout  << "operator== start"<<std::endl;
   if (this->variance() != other.variance()) throw std::logic_error("mismatching co/contravariance");
   if (this->m_size != other.m_size) throw std::logic_error("mismatching lengths");
   constexpr bool pr=false;//bool pr=m_size<3;
//   bool pr=m_mpi_rank>0;
   int diff=0;
   if (!m_cache.io && !other.m_cache.io) { // std::cout << "both in memory"<<std::endl;
    size_t off = (m_replicated && !other.m_replicated) ? other.m_segment_offset : 0;
    size_t otheroff = (!m_replicated && other.m_replicated) ? m_segment_offset : 0;
     for (size_t i=0; i<std::min(m_segment_length,other.m_segment_length); i++) {
      diff = diff || m_cache.buffer[off+i] != other.m_cache.buffer[otheroff+i];
           if (pr) std::cout <<m_mpi_rank<<" other.buffer["<<otheroff+i<<"]="<<other.m_cache.buffer[otheroff+i]<<std::endl;
           if (pr) std::cout <<m_mpi_rank<<" buffer["<<off+i<<"]="<<m_cache.buffer[off+i]<<std::endl;
     }
   } else if (m_cache.io && other.m_cache.io) { // std::cout << "both cached"<<std::endl;
    size_t cachelength = std::min(m_cache.preferred_length,other.m_cache.preferred_length);
//    std::cout << m_mpi_rank<<"cachelength="<<cachelength<<std::endl;
    size_t off=0, otheroff=0;
    //     m_cache.ensure( (m_replicated == other.m_replicated) ? 0: other.m_segment_offset); std::cout << "m_cache.length="<<m_cache.length<<std::endl;
    //     other.m_cache.ensure( (m_replicated == other.m_replicated) ? 0: m_segment_offset); std::cout << "other.m_cache.length="<<other.m_cache.length<<std::endl;
    if (pr) std::cout << "move cache to "<< ((!m_replicated ||  other.m_replicated) ? 0: other.m_segment_offset)<<std::endl;
    if (pr) std::cout << "move other cache to "<<((m_replicated || !other.m_replicated) ? 0: m_segment_offset)<<std::endl;
    m_cache.move( (!m_replicated ||  other.m_replicated) ? 0: other.m_segment_offset, cachelength);
    other.m_cache.move((m_replicated || !other.m_replicated) ? 0: m_segment_offset, cachelength);
    while(m_cache.length && other.m_cache.length && off < m_segment_length && otheroff < other.m_segment_length ) {
     if (pr) std::cout <<m_mpi_rank<< " buffer to copy in range "<<m_cache.offset<<"="<<other.m_cache.offset<<" for "<<m_cache.length<<"="<<other.m_cache.length<<std::endl;
     if (pr) for (size_t i=0; i<m_cache.length; i++) std::cout <<" "<<other.m_cache.buffer[otheroff+i]; std::cout <<std::endl;
     for (size_t i=0; i<m_cache.length; i++) {
      diff = diff || m_cache.buffer[off+i] != other.m_cache.buffer[otheroff+i];
           if (pr) std::cout <<m_mpi_rank<<" buffer["<<off+i<<"]="<<m_cache.buffer[off+i]<<std::endl;
     }
     m_cache.dirty=true;
     if (m_cache.io) ++m_cache; else off+=cachelength;
     if (other.m_cache.io) ++other.m_cache; else otheroff+=cachelength;
     if (pr) std::cout <<m_mpi_rank<< "end of cache loop, off="<<off<<", otheroff="<<otheroff<<std::endl;
    }
   } else if (m_cache.io) { // std::cout<< "source in memory, result cached"<<std::endl;
    size_t cachelength = m_cache.preferred_length;
    size_t off=0;
    size_t otheroff=(m_replicated == other.m_replicated) ? 0: m_segment_offset;
    m_cache.move( (m_replicated == other.m_replicated) ? 0: other.m_segment_offset, cachelength );
    other.m_cache.move( 0, other.m_segment_length);
    while(m_cache.length && off < m_segment_length && otheroff < other.m_segment_length ) {
//      std::cout <<m_mpi_rank<< " buffer to copy in range "<<m_cache.offset<<"="<<other.m_cache.offset<<" for "<<m_cache.length<<"="<<other.m_cache.length<<std::endl;
//      for (size_t i=0; i<m_cache.length; i++) std::cout <<" "<<other.m_cache.buffer[otheroff+i]; std::cout <<std::endl;
//      std::cout <<"off = "<<off<<std::endl;
     for (size_t i=0; i<m_cache.length; i++) {
      diff = diff || m_cache.buffer[off+i] != other.m_cache.buffer[otheroff+i];
      //     if (pr) std::cout <<m_mpi_rank<<" buffer["<<off+i<<"]="<<m_cache.buffer[off+i]<<std::endl;
     }
      m_cache.dirty=true;
     ++m_cache;
     otheroff+=cachelength;
     if (pr) std::cout << "end of cache loop, off="<<off<<", otheroff="<<otheroff<<std::endl;
    }
   } else if (other.m_cache.io) {  // std::cout << "source cached, result in memory"<<std::endl;
    size_t cachelength = other.m_cache.preferred_length;
    size_t otheroff=0;
    size_t off=(m_replicated == other.m_replicated) ? 0: other.m_segment_offset;
    m_cache.move(0, m_segment_length);
    other.m_cache.move((m_replicated == other.m_replicated) ? 0: m_segment_offset, cachelength );
    while(other.m_cache.length && off < m_segment_length && otheroff < other.m_segment_length ) {
     if (pr) std::cout <<m_mpi_rank<< " buffer to copy in range "<<m_cache.offset<<"="<<other.m_cache.offset<<" for "<<m_cache.length<<"="<<other.m_cache.length<<std::endl;
     //    for (size_t i=0; i<m_cache.length; i++) std::cout <<" "<<other.m_cache.buffer[otheroff+i]; std::cout <<std::endl;
     for (size_t i=0; i<other.m_cache.length; i++) {
      diff = diff || m_cache.buffer[off+i] != other.m_cache.buffer[otheroff+i];
      //     if (pr) std::cout <<m_mpi_rank<<" buffer["<<off+i<<"]="<<m_cache.buffer[off+i]<<std::endl;
     }
     off+=cachelength;
     ++other.m_cache;
     if (pr) std::cout << "end of cache loop, off="<<off<<", otheroff="<<otheroff<<std::endl;
    }
   }
//   if (pr) std::cout << "operator= after copy loop, other="<<other<<std::endl;
//   if (pr) std::cout << "operator= after copy loop, this="<<*this<<std::endl;
#ifdef USE_MPI
   if (!m_replicated && ! other.m_replicated) {
    for (int rank=0; rank < m_mpi_size; rank++) {
      MPI_Allreduce(&diff,&diff,1,MPI_INT,MPI_SUM,m_communicator);
    }
   }
#endif
//std::cout  << "operator== end"<<std::endl;
   return diff==0;
  }

 };
 template <class scalar>
 inline std::ostream& operator<<(std::ostream& os, PagedVector<scalar> const& obj) { return os << obj.str(); }

}
#endif // PAGEDVECTOR_H
