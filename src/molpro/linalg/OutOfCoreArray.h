#ifndef OUTOFCOREARRAY_H
#define OUTOFCOREARRAY_H
#ifdef TIMING
#include <chrono>
#endif
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <climits>
#include <cstring>
#include <cstddef>
#include <cstdlib>
#include <unistd.h>
#include <cassert>
#include <fstream>
#include <sstream>
#include <ostream>
#include <map>
#include <vector>
#include <iostream>


#ifdef HAVE_MPI_H
#include <mpi.h>
#ifdef HAVE_PPIDD_H
#include "ppidd.h"
#endif
#ifndef MPI_COMM_COMPUTE
#ifdef HAVE_PPIDD_H
#define MPI_COMM_COMPUTE MPI_Comm_f2c(PPIDD_Worker_comm())
#else
#define MPI_COMM_COMPUTE MPI_COMM_WORLD
#endif
#endif
#else
#define MPI_COMM_COMPUTE 0
using MPI_Comm = int;
#endif

namespace molpro {
namespace linalg {

/*!
  * \brief A class that implements a fixed-size vector container that has the following public features:
  * - data held on backing store, distributed over MPI ranks, instead of in memory
  * - alternatively, mapping of an externally-owned buffer
  * - opaque implementation of BLAS including dot(), axpy(), scal()
  * - additional BLAS overloads to combine with a simple sparse vector represented in std::map
  * - efficient import and export of data ranges
  * - additional functions to find the largest elements, and to print the vector
  * - underlying type of elements
  * \tparam T the type of elements of the vector
  * \tparam Allocator alternative to std::allocator
  */
template<class T=double,
        size_t default_buffer_size = 102400,
        class Allocator =
#ifdef MEMORY_MEMORY_H
        memory::allocator<T>
#else
        std::allocator<T>
#endif
        >
class OutOfCoreArray {
  typedef double scalar_type; //TODO implement this properly from T
 public:
  typedef T value_type;
  /*!
   * \brief Construct an object without any data.
   */
  explicit OutOfCoreArray(size_t length = 0, MPI_Comm mpi_communicator = MPI_COMM_COMPUTE)
      : m_size(length),
        m_communicator(mpi_communicator), m_mpi_size(mpi_size()), m_mpi_rank(mpi_rank()),
        m_data(nullptr),
        m_sync(true),
        m_segment_offset(seg_offset()),
        m_segment_length(seg_length()),
        m_buf_size(0),
        m_file(make_file()) {
  }
  /*!
   * @brief Copy constructor
   * @param source
   * @param option (to comply with other vector classes)
   */
  OutOfCoreArray(const OutOfCoreArray& source, int option=0)
      : m_size(source.m_size),
        m_communicator(source.m_communicator), m_mpi_size(mpi_size()), m_mpi_rank(mpi_rank()),
        m_data(nullptr),
        m_sync(source.m_sync),
        m_segment_offset(seg_offset()),
        m_segment_length(seg_length()),
        m_buf_size(m_segment_length > default_buffer_size ? default_buffer_size : m_segment_length),
        m_file(make_file()) {
#ifdef TIMING
    auto start=std::chrono::steady_clock::now();
#endif
    if (source.m_data == nullptr) {
        throw std::logic_error("Source object must be pointing to external buffer");
    }
    if (!source.m_sync
        and option == 0 // this is a fix to condone IterativeSolver playing tricks where synchronisation isn't actually needed
        )
      throw std::logic_error("Source object must be synchronised across ranks"); //or just sync this?
    if (m_buf_size == 0) return;
    size_t offset = 0;
    while ((m_segment_length-offset)/m_buf_size != 0) {
        m_file.seekp(offset * sizeof(T));
        m_file.write((const char*) (source.m_data + m_segment_offset + offset), m_buf_size * sizeof(T));
        offset += m_buf_size;
    }
    if (offset < m_segment_length) {
        size_t rest = m_segment_length-offset;
        m_file.seekp(offset * sizeof(T));
        m_file.write((const char*) (source.m_data + m_segment_offset + offset), rest * sizeof(T));
    }
    // close file?? re-wind??
    // makes put() obsolete?
#ifdef TIMING
    auto end=std::chrono::steady_clock::now();
    std::cout <<" OutOfCoreArray copy constructor length="<<m_size
    <<", seconds="<<std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count()*1e-9
        <<", bandwidth="<<m_size*8*1e3/std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count()<<"MB/s"
    <<std::endl;
#endif
  }

  /*!
  * @brief Construct an object mapped from external data
   *
   * @param buffer
   * @param length
   * @param mpi_communicator Not used by this object, but can be specified so that it is taken when copy-constructing from this object
   */

  OutOfCoreArray(T* buffer, size_t length, MPI_Comm mpi_communicator = MPI_COMM_COMPUTE)
      : m_size(length),
        m_communicator(mpi_communicator), m_mpi_size(mpi_size()), m_mpi_rank(mpi_rank()),
        m_segment_offset(seg_offset()),
        m_segment_length(seg_length()),
        m_data(buffer == nullptr ? nullptr : buffer),
        m_sync(true),
        m_buf_size(0),
        m_file(std::fstream()) {
      if (m_data == nullptr) {
          throw std::logic_error("Source object must be pointing to external buffer");
      }
//    std::cout << "OutOfCoreArray map constructor "<<buffer<<" : "<<&(*this)[0]<<std::endl;
  }

  //~OutOfCoreArray() = default;
  ~OutOfCoreArray() {
      m_file.close();
  }

  bool synchronised() const { return m_sync; }

 private:
  std::fstream make_file() {
    std::fstream file;
    char* tmpname = strdup("tmpfileXXXXXX");
    {
      auto ifile =
          mkstemp(tmpname); // actually create the file atomically with the name to avoid race resulting in duplicates or conflicts
      if (ifile < 0) throw std::runtime_error(std::string("Cannot open cache file ") + tmpname);
      close(ifile);
    }
    file.open(tmpname, std::ios::out | std::ios::in | std::ios::binary);
    if (!file.is_open()) throw std::runtime_error(std::string("Cannot open cache file ") + tmpname);
    unlink(tmpname);
    free(tmpname);
    return file;
  }

  int mpi_size() {
#ifdef HAVE_MPI_H
    int result;
    MPI_Initialized(&result);
#ifdef HAVE_PPIDD_H
    if (!result) PPIDD_Initialize(0, nullptr, PPIDD_IMPL_DEFAULT);
#else
    if (!result) MPI_Init(0,nullptr);
#endif
    MPI_Comm_size(m_communicator, &result);
#else
    int result = 1;
#endif
    return result;
  }
  int mpi_rank() {
#ifdef HAVE_MPI_H
    int result;
    MPI_Comm_rank(m_communicator, &result);
#else
    int result = 0;
#endif
    return result;
  }

  size_t m_size; //!< How much data
  const MPI_Comm m_communicator; //!< the MPI communicator for distributed data
  int m_mpi_size;
  int m_mpi_rank;
  T* m_data; //!< the data buffer (if in an externally-provided buffer) or nullptr (if owned by this object)
  bool m_sync; //!< whether replicated data is in sync on all processes
  size_t m_segment_offset; //!< offset in the overall data object of this process' data
  size_t m_segment_length; //!< length of this process' data
  size_t m_buf_size; //!< buffer size
  mutable std::fstream m_file;

  /*!
 * Calculate the segment length for asynchronous operations on replicated vectors
 * @return the length of a vector segment
 */
  size_t seg_length() const {
    size_t alpha = m_size / m_mpi_size;
    size_t beta = m_size % m_mpi_size;
    size_t segl = (m_mpi_rank < beta) ? alpha + 1 : alpha;
    return segl;
  }

  /*!
 * Calculate the segment offset for asynchronous operations on replicated vectors
 * @return the offset of a vector segment
 */
  size_t seg_offset() const {
    size_t alpha = m_size / m_mpi_size;
    size_t beta = m_size % m_mpi_size;
    size_t segoff = (m_mpi_rank < beta) ? (alpha + 1) * m_mpi_rank : alpha * m_mpi_rank + beta;
    return segoff;
  }

  /*!
 * Get a pointer to the entire data structure if possible
 * @return the pointer, or nullptr if caching/distribution means that the whole vector is not available
 */
  T* data() {
    return m_data;
  }

 public:
  /*!
   * \brief Update a range of the object data with the contents of a provided buffer
   * \param buffer The provided data segment
   * \param length Length of data
   * \param offset Offset in the object of buffer[0]
   */
  void put(const T* buffer, size_t length, size_t offset) {
    if (m_data != nullptr) { // data in buffer
      std::copy(buffer, buffer + length, m_data + offset);
      return;
    }
    if (offset > m_segment_offset + m_segment_length or offset + length < m_segment_offset)
      return;
    size_t buffer_offset = 0;
    if (offset < m_segment_offset) {
      buffer_offset = m_segment_offset - offset;
      offset = m_segment_offset;
      length -= m_segment_offset - offset;
    }
    if (offset + length > std::min(m_size, m_segment_offset + m_segment_length))
      length = std::min(m_size,
                        m_segment_offset + m_segment_length) - offset;
    m_file.seekp(offset * sizeof(T));
    m_file.write((const char*) buffer + buffer_offset, length * sizeof(T));
  }

  /*!
   * \brief Read a range of the object data into a provided buffer
   * @param buffer
   * @param length
   * @param offset
   */
  void get(T* buffer, size_t length, size_t offset) const {
    if (m_data != nullptr) { // data in buffer
      std::copy(m_data + offset, m_data + offset + length, buffer);
      return;
    } //otherwise - data distributed so use local offset
    if (offset >= m_segment_length or offset + length < 0)
      return; // through exception??
    size_t buffer_offset = 0;
    if (offset < 0) {
      buffer_offset = -offset;
      offset = 0;
      length -= buffer_offset;
    }
    if ((offset + length) > m_segment_length)
      length = m_segment_length - offset;
    m_file.seekg(offset * sizeof(T));
    m_file.read((char*) buffer + buffer_offset, length * sizeof(T));
  }

    /*!
     * @brief Return a reference to an element of the data
     * @param pos Offset of the data
     * @return
     */
  const T operator[](size_t pos) const {
      if (m_data != nullptr) {
          return m_data[pos];
      } else {
          T value;
          int root_rank;
          size_t displs[m_mpi_size];
          for (int i = 0; i < m_mpi_size; i++) {
              displs[i] = (i < m_size%m_mpi_size) ? (m_size/m_mpi_size+1)*i : (m_size/m_mpi_size)*i + m_size%m_mpi_size;
              if (pos < displs[i]) {
                 root_rank = i-1;
                 break;
              } else if (i == (m_mpi_size - 1)) {
                 root_rank = i;
              }
          }
          if (root_rank == m_mpi_rank) // only the owner rank
             get(&value,1,pos - m_segment_offset);
#ifdef HAVE_MPI_H
          MPI_Bcast(&value, 1, MPI_DOUBLE, root_rank, m_communicator);
#endif
          return value;
      }
  }

  /*!
   * @brief Return the number of elements of data
   * @return
   */
  size_t size() const { return m_size; }

  /*!
   * @brief Produce a printable representation of the object
   * @param verbosity How much to print
   * @param columns Width in characters
   * @return
   */
  std::string str(int verbosity = 0, unsigned int columns = UINT_MAX) const {
    std::ostringstream os;
    os << "OutOfCoreArray object:";
    if (m_data != nullptr) {
        for (size_t k = 0; k < m_segment_length; k++)
            os << " " << (*this)[m_segment_offset + k];
    } else if (m_buf_size > 0) {
        std::vector<T> buffer;
        buffer.reserve(m_buf_size);
        size_t offset = 0;
        while ((m_segment_length-offset)/m_buf_size != 0) {
            get(&buffer[0],m_buf_size,offset);
            for (size_t i = 0; i < m_buf_size; i++)
                os << " " << buffer[i];
            offset += m_buf_size;
        }
        if (offset < m_segment_length) {
            size_t rest = m_segment_length-offset;
            get(&buffer[0],rest,offset);
            for (size_t i = 0; i < rest; i++)
                os << " " << buffer[i];
        }
    }
    os << std::endl;
    return os.str();
  }

  /*
   * \brief replicate the vectors over the MPI ranks
   * \return
   */
  void sync() {
#ifdef HAVE_MPI_H
    //std::cout <<m_mpi_rank<<"before broadcast this="<<*this<<std::endl;
      if (m_data != nullptr && !m_sync) { //only case when can be out of sync
          size_t alpha = m_size/m_mpi_size;
          size_t beta = m_size%m_mpi_size;
          int chunks[m_mpi_size];
          int displs[m_mpi_size];
          for (int i = 0; i < m_mpi_size; i++) {
              chunks[i] = (i < beta) ? alpha+1 : alpha;
              displs[i] = (i < beta) ? chunks[i]*i : chunks[i]*i + beta;
              //            if (m_mpi_rank == 0) std::cout<<"Chunk: "<<chunks[i]<<" Displ: "<<displs[i]<<std::endl;
          }
          MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,&m_data[0],chunks,displs,MPI_DOUBLE,m_communicator);
      }
#endif
      m_sync = true;
  }

  /*!
   * \brief Add a constant times another object to this object
   * \param a The factor to multiply.
   * \param other The object to be added to this.
   * \return
   */
  void axpy(scalar_type a, const OutOfCoreArray<T,default_buffer_size>& other) {
    if (m_size != other.m_size) throw std::logic_error("mismatching lengths"); //m_size -> other.m_size - correct?
    if (m_data != nullptr) {
        if (other.m_data == nullptr && other.m_buf_size > 0) {
            std::vector<T> buffer;
            buffer.reserve(other.m_buf_size);
            size_t offset = 0;
            while ((m_segment_length-offset)/other.m_buf_size != 0) {
                other.get(&buffer[0],other.m_buf_size,offset);
                for (size_t i = 0; i < other.m_buf_size; i++)
                    m_data[m_segment_offset + offset + i] += a * buffer[i]; //or use BLAS
                offset += other.m_buf_size;
            }
            if (offset < m_segment_length) {
                size_t rest = m_segment_length-offset;
                other.get(&buffer[0],rest,offset);
                for (size_t i = 0; i < rest; i++)
                    m_data[m_segment_offset + offset + i] += a * buffer[i]; //or use BLAS
            } // Too lengthy? Wrap in a function?
        } else { // Not expected!? Should check if the same vector => "==" operator needed??
            for (size_t i = 0; i < m_segment_length; i++) // check that of same length?
                m_data[m_segment_offset + i] += a * other.m_data[other.m_segment_offset + i]; //or use BLAS
        }
        m_sync = false;
    } else { // Not expected!?
        throw std::logic_error("Axpy should not be called for a vector written to disk");
    }
  }

  /*!
    * \brief Add a constant times a sparse vector to this object
    * \param a The factor to multiply.
    * \param other The object to be added to this.
    * \return
    */
  void axpy(scalar_type a, const std::map<size_t, T>& other) {
    if (m_data == nullptr) {
        throw std::logic_error("Axpy should not be called for a vector written to disk");
    }
    for (const auto& o: other)
      if (o.first >= m_segment_offset && o.first < m_segment_offset + m_segment_length)
        m_data[o.first] += a * o.second;
    m_sync = false;
  }

  /*!
   * \brief Scalar product of two objects.
   * \param other The object to be contracted with this.
   * \return
   */
  scalar_type dot(const OutOfCoreArray<T,default_buffer_size>& other) const {
    if (this->m_size != other.m_size) throw std::logic_error("mismatching lengths");
    scalar_type result = 0;
    if (this == &other) {
        if (m_data != nullptr) {
            for (size_t i = 0; i < m_segment_length; i++)
                result += m_data[m_segment_offset + i] * m_data[m_segment_offset + i]; //or use BLAS
        } else if (m_buf_size > 0) {
            std::vector<T> buffer;
            buffer.reserve(m_buf_size);
            size_t offset = 0;
            while ((m_segment_length-offset)/m_buf_size != 0) {
                get(&buffer[0],m_buf_size,offset);
                for (size_t i = 0; i < m_buf_size; i++)
                    result += buffer[i] * buffer[i]; //or use BLAS
                offset += m_buf_size;
            }
            if (offset < m_segment_length) {
                size_t rest = m_segment_length-offset;
                get(&buffer[0],rest,offset);
                for (size_t i = 0; i < rest; i++)
                    result += buffer[i] * buffer[i]; //or use BLAS
            }
        }
    } else {
        if (m_data != nullptr) {
            if (other.m_data != nullptr) {
                for (size_t i = 0; i < m_segment_length; i++) // check m_segment_length == other.m_segment_length?
                    result += m_data[m_segment_offset + i] * other.m_data[other.m_segment_offset + i]; //or use BLAS
            } else if (other.m_buf_size > 0) {
                std::vector<T> buffer;
                buffer.reserve(other.m_buf_size);
                size_t offset = 0;
                while ((m_segment_length-offset)/other.m_buf_size != 0) {
                    other.get(&buffer[0],other.m_buf_size,offset);
                    for (size_t i = 0; i < other.m_buf_size; i++)
                        result += m_data[m_segment_offset + offset + i] * buffer[i]; //or use BLAS
                    offset += other.m_buf_size;
                }
                if (offset < m_segment_length) {
                    size_t rest = m_segment_length-offset;
                    other.get(&buffer[0],rest,offset);
                    for (size_t i = 0; i < rest; i++)
                        result += m_data[m_segment_offset + offset + i] * buffer[i]; //or use BLAS
                }
            }
        } else if (m_buf_size > 0) {
            if (other.m_data != nullptr) {
                std::vector<T> buffer;
                buffer.reserve(m_buf_size);
                size_t offset = 0;
                while ((m_segment_length-offset)/m_buf_size != 0) {
                    get(&buffer[0],m_buf_size,offset);
                    for (size_t i = 0; i < m_buf_size; i++)
                        result += buffer[i] * other.m_data[other.m_segment_offset + offset + i]; //or use BLAS
                    offset += m_buf_size;
                }
                if (offset < m_segment_length) {
                    size_t rest = m_segment_length-offset;
                    get(&buffer[0],rest,offset);
                    for (size_t i = 0; i < rest; i++)
                        result += buffer[i] * other.m_data[other.m_segment_offset + offset + i]; //or use BLAS
                }
            } else {
                std::vector<T> lbuffer,rbuffer;
                lbuffer.reserve(m_buf_size);
                rbuffer.reserve(other.m_buf_size); //check same length??
                size_t offset = 0;
                while ((m_segment_length-offset)/m_buf_size != 0) {
                    get(&lbuffer[0],m_buf_size,offset);
                    other.get(&rbuffer[0],other.m_buf_size,offset);
                    for (size_t i = 0; i < m_buf_size; i++)
                        result += lbuffer[i] * rbuffer[i]; //or use BLAS
                    offset += m_buf_size;
                }
                if (offset < m_segment_length) {
                    size_t rest = m_segment_length-offset;
                    get(&lbuffer[0],rest,offset);
                    other.get(&rbuffer[0],rest,offset);
                    for (size_t i = 0; i < rest; i++)
                        result += lbuffer[i] * rbuffer[i]; //or use BLAS
                }
            }
        }
    }
#ifdef HAVE_MPI_H
    //std::cout <<m_mpi_rank<<" dot result before reduce="<<result<<std::endl;
    double resultLocal = result;
    double resultGlobal = result;
    MPI_Allreduce(&resultLocal, &resultGlobal, 1, MPI_DOUBLE, MPI_SUM, m_communicator);
    result = resultGlobal;
    // std::cout <<m_mpi_rank<<" dot result after reduce="<<result<<std::endl;
#endif
    return result;
  }

  /*!
   * \brief Scalar product with a sparse vector
   * \param other The object to be contracted with this.
   * \return
   */
  scalar_type dot(const std::map<size_t, T>& other) const {
    scalar_type result = 0;
    if (m_data != 0) {
        for (const auto& o: other)
            if (o.first >= m_segment_offset && o.first < m_segment_offset + m_segment_length)
                result += o.second * m_data[o.first];
    } else if (m_buf_size > 0) {
        std::vector<T> buffer;
        buffer.reserve(m_buf_size);
        size_t offset = 0;
        while ((m_segment_length-offset)/m_buf_size != 0) {
            get(&buffer[0],m_buf_size,offset);
            for (const auto& o: other)
                if ((o.first >= m_segment_offset + offset) && (o.first < m_segment_offset + offset + m_buf_size))
                    result += o.second * buffer[o.first - m_segment_offset - offset]; //or use BLAS
            offset += m_buf_size;
        }
        if (offset < m_segment_length) {
            size_t rest = m_segment_length-offset;
            get(&buffer[0],rest,offset);
            for (const auto& o: other)
                if ((o.first >= m_segment_offset + offset) && (o.first < m_segment_offset + offset + rest))
                    result += o.second * buffer[o.first - m_segment_offset - offset]; //or use BLAS
        }
    }
#ifdef HAVE_MPI_H
    double resultLocal = result;
    double resultGlobal = result;
    MPI_Allreduce(&resultLocal, &resultGlobal, 1, MPI_DOUBLE, MPI_SUM, m_communicator); // FIXME needs attention for non-double
    result = resultGlobal;
#endif
    return result;
  }

  /*!
     * \brief scal Scale the object by a factor.
     * \param a The factor to scale by. If a is zero, then the current contents of the object are ignored.
     */
  void scal(scalar_type a) {
      if (m_data == nullptr) {
          throw std::logic_error("Scal() should not be called for a vector written to disk");
      }
      if (a != 0) {
          for (size_t i = 0; i < m_segment_length; i++)
              m_data[m_segment_offset + i] *= a;
      } else {
          for (size_t i = 0; i < m_segment_length; i++)
              m_data[m_segment_offset + i] = 0;
      }
      m_sync = false;
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
  std::tuple<std::vector<size_t>, std::vector<T> > select(
      const OutOfCoreArray<T,default_buffer_size>& measure,
      const size_t maximumNumber = 1000,
      const scalar_type threshold = 0
  ) const {
    typedef std::multimap<T, size_t, std::greater<T> > MyMap;
    MyMap sortlist;
    const auto& measur = dynamic_cast <const OutOfCoreArray<T,default_buffer_size>&> (measure); //Correct??
    if (m_size != measur.m_size) throw std::logic_error("mismatching lengths");
    if (m_data == nullptr && measur.m_data != nullptr) throw std::logic_error("mismatching distribution status");
    if (m_data != nullptr && measur.m_data == nullptr) throw std::logic_error("mismatching distribution status");
    if (m_data != nullptr) {
        if (this == &measur) {
            for (size_t i = 0; i < m_segment_length; i++) {
                auto test = m_data[m_segment_offset + i] * m_data[m_segment_offset + i];
                if (test < 0) test = -test; //?
                if (test > threshold) {
                    sortlist.insert(std::make_pair(test, i + m_segment_offset));
                    if (sortlist.size() > maximumNumber) sortlist.erase(std::prev(sortlist.end()));
                }
            }
        } else {
            for (size_t i = 0; i < m_segment_length; i++) {
                auto test = m_data[m_segment_offset + i] * measur.m_data[measur.m_segment_offset + i];
                if (test < 0) test = -test; //?
                if (test > threshold) {
                    sortlist.insert(std::make_pair(test, i + m_segment_offset));
                    if (sortlist.size() > maximumNumber) sortlist.erase(std::prev(sortlist.end()));
                }
            }
        }
    } else if (m_buf_size > 0) {
        //throw std::logic_error("vectors on disk not expected in select()"); //correct?
        if (this == &measur) {
            std::vector<T> buffer;
            buffer.reserve(m_buf_size);
            size_t offset = 0;
            while ((m_segment_length-offset)/m_buf_size != 0) {
                get(&buffer[0],m_buf_size,offset);
                for (size_t i = 0; i < m_buf_size; i++) {
                    auto test = buffer[i] * buffer[i]; //or use BLAS
                    if (test < 0) test = -test;
                    if (test > threshold) {
                        sortlist.insert(std::make_pair(test, i + m_segment_offset + offset));
                        if (sortlist.size() > maximumNumber) sortlist.erase(std::prev(sortlist.end()));
                    }
                }
                offset += m_buf_size;
            }
            if (offset < m_segment_length) {
                size_t rest = m_segment_length-offset;
                get(&buffer[0],rest,offset);
                for (size_t i = 0; i < rest; i++) {
                    auto test = buffer[i] * buffer[i]; //or use BLAS
                    if (test < 0) test = -test;
                    if (test > threshold) {
                        sortlist.insert(std::make_pair(test, i + m_segment_offset + offset));
                        if (sortlist.size() > maximumNumber) sortlist.erase(std::prev(sortlist.end()));
                    }
                }
            }
        } else {
            std::vector<T> lbuffer,rbuffer;
            lbuffer.reserve(m_buf_size);
            rbuffer.reserve(measur.m_buf_size);
            size_t offset = 0;
            while ((m_segment_length-offset)/m_buf_size != 0) {
                get(&lbuffer[0],m_buf_size,offset);
                measur.get(&rbuffer[0],measur.m_buf_size,offset);
                for (size_t i = 0; i < m_buf_size; i++) {
                    auto test = lbuffer[i] * rbuffer[i]; //or use BLAS
                    if (test < 0) test = -test;
                    if (test > threshold) {
                        sortlist.insert(std::make_pair(test, i + m_segment_offset + offset));
                        if (sortlist.size() > maximumNumber) sortlist.erase(std::prev(sortlist.end()));
                    }
                }
                offset += m_buf_size;
            }
            if (offset < m_segment_length) {
                size_t rest = m_segment_length-offset;
                get(&lbuffer[0],rest,offset);
                measur.get(&rbuffer[0],rest,offset);
                for (size_t i = 0; i < rest; i++) {
                    auto test = lbuffer[i] * rbuffer[i]; //or use BLAS
                    if (test < 0) test = -test;
                    if (test > threshold) {
                        sortlist.insert(std::make_pair(test, i + m_segment_offset + offset));
                        if (sortlist.size() > maximumNumber) sortlist.erase(std::prev(sortlist.end()));
                    }
                }
            }
        }
    }
    while (sortlist.size() > maximumNumber) sortlist.erase(std::prev(sortlist.end())); //redundant?
    std::vector<size_t> indices;
    indices.reserve(sortlist.size());
    std::vector<T> values;
    values.reserve(sortlist.size());
    std::transform(sortlist.begin(), sortlist.end(), std::back_inserter(indices), //or indices.begin()?
                   [](const typename MyMap::value_type& val){return val.second;} );
    std::transform(sortlist.begin(), sortlist.end(), std::back_inserter(values),
                   [](const typename MyMap::value_type& val){return val.first;} );
#ifdef HAVE_MPI_H
    std::vector<int> intindices(indices.begin(),indices.end());
    std::vector<int> long_intindices;
    long_intindices.reserve(maximumNumber * m_mpi_size);
    std::vector<T> long_values;
    long_values.reserve(maximumNumber * m_mpi_size);
    int chunks[m_mpi_size],displs[m_mpi_size];
    int ssize = sortlist.size();
    MPI_Allgather(&ssize,1,MPI_INT,chunks,1,MPI_INT,m_communicator);
    displs[0] = 0;
    for (int i = 1; i < m_mpi_size; i++) {
        displs[i] = chunks[i-1];
    }
    MPI_Allgatherv(values.data(), values.size(), MPI_DOUBLE,
            long_values.data(), chunks, displs, MPI_DOUBLE, m_communicator); //non-blocking??
    MPI_Allgatherv(intindices.data(), intindices.size(), MPI_INT,
                  long_intindices.data(), chunks, displs, MPI_INT, m_communicator); //non-blocking??
    std::vector<size_t> long_indices(long_intindices.begin(),long_intindices.end());
    MyMap fsortlist;
    for (size_t i = 0; i < long_indices.size(); i++) {
        fsortlist.insert(std::make_pair(long_values[i], long_indices[i]));
        if (fsortlist.size() > maximumNumber) fsortlist.erase(std::prev(fsortlist.end()));
    }
    while (fsortlist.size() > maximumNumber) fsortlist.erase(std::prev(fsortlist.end())); //redundant?
    std::transform(fsortlist.begin(), fsortlist.end(), indices.begin(),
                   [](const typename MyMap::value_type& val){return val.second;} );
    std::transform(fsortlist.begin(), fsortlist.end(), values.begin(),
                   [](const typename MyMap::value_type& val){return val.first;} );
#endif
    return std::make_tuple(indices, values);
  };

  /*!
   * \brief Copy from one object to another, adjusting size if needed.
   * \param other The source of data.
   * \return
   */
  OutOfCoreArray& operator=(const OutOfCoreArray& other) {
      if (m_size != other.m_size) throw std::logic_error("Two objects must have same size");
      if (m_data == nullptr) throw std::logic_error("Can't change vector on disk");
      if (other.m_data != nullptr) {
          for (size_t i = 0; i < m_segment_length; i++)
              m_data[i + m_segment_offset] = other.m_data[i + other.m_segment_offset];
      } else if (other.m_buf_size > 0) {
          std::vector<T> buffer;
          buffer.reserve(other.m_buf_size);
          size_t offset = 0;
          while ((other.m_segment_length-offset)/other.m_buf_size != 0) {
              other.get(&buffer[0],other.m_buf_size,offset);
              for (size_t i = 0; i < other.m_buf_size; i++)
                  m_data[m_segment_offset + offset + i] = buffer[i]; //or use BLAS
              offset += other.m_buf_size;
          }
          if (offset < other.m_segment_length) {
              size_t rest = other.m_segment_length-offset;
              other.get(&buffer[0],rest,offset);
              for (size_t i = 0; i < rest; i++)
                  m_data[m_segment_offset + offset + i] = buffer[i]; //or use BLAS
          }
      }
      return *this;
  }

  /*!
   * @brief Test for equality of two objects
   * @param other
   * @return
   */
  bool operator==(const OutOfCoreArray& other) const {
    if (m_size != other.m_size) throw std::logic_error("mismatching lengths");
    int diff = 0;
    if (m_data != nullptr && other.m_data != nullptr) { //both in ext. buffer
        for (size_t i = 0; i < m_segment_length; i++)
            diff = diff || m_data[i + m_segment_offset] != other.m_data[i + other.m_segment_offset];
    } else if (m_data == nullptr && other.m_data == nullptr && m_buf_size > 0) { // both dist., on disk
        std::vector<T> lbuffer,rbuffer;
        lbuffer.reserve(m_buf_size);
        rbuffer.reserve(other.m_buf_size);
        size_t offset = 0;
        while ((m_segment_length-offset)/m_buf_size != 0) {
            get(&lbuffer[0],m_buf_size,offset);
            other.get(&rbuffer[0],other.m_buf_size,offset);
            for (size_t i = 0; i < m_buf_size; i++)
                diff = diff || lbuffer[i] != rbuffer[i]; //or use BLAS
            offset += m_buf_size;
        }
        if (offset < m_segment_length) {
            size_t rest = m_segment_length-offset;
            get(&lbuffer[0],rest,offset);
            other.get(&rbuffer[0],rest,offset);
            for (size_t i = 0; i < rest; i++)
                diff = diff || lbuffer[i] != rbuffer[i]; //or use BLAS
        }
    } else if (m_data != nullptr && other.m_buf_size > 0) { // this in ext. buffer, other dist., on disk
        std::vector<T> buffer;
        buffer.reserve(other.m_buf_size);
        size_t offset = 0;
        while ((m_segment_length-offset)/other.m_buf_size != 0) {
            other.get(&buffer[0],other.m_buf_size,offset);
            for (size_t i = 0; i < other.m_buf_size; i++)
                diff = diff || m_data[m_segment_offset + offset + i] != buffer[i]; //or use BLAS
            offset += other.m_buf_size;
        }
        if (offset < m_segment_length) {
            size_t rest = m_segment_length-offset;
            other.get(&buffer[0],rest,offset);
            for (size_t i = 0; i < rest; i++)
                diff = diff || m_data[m_segment_offset + offset + i] != buffer[i]; //or use BLAS
        }
    } else if (m_buf_size > 0) { // this dist., on disk, other in ext. buffer
        std::vector<T> buffer;
        buffer.reserve(m_buf_size);
        size_t offset = 0;
        while ((m_segment_length-offset)/m_buf_size != 0) {
            get(&buffer[0],m_buf_size,offset);
            for (size_t i = 0; i < m_buf_size; i++)
                diff = diff || buffer[i] != other.m_data[other.m_segment_offset + offset + i]; //or use BLAS
            offset += m_buf_size;
        }
        if (offset < m_segment_length) {
            size_t rest = m_segment_length-offset;
            get(&buffer[0],rest,offset);
            for (size_t i = 0; i < rest; i++)
                diff = diff || buffer[i] != other.m_data[other.m_segment_offset + offset + i]; //or use BLAS
        }
    }
#ifdef HAVE_MPI_H
    MPI_Allreduce(MPI_IN_PLACE, &diff, 1, MPI_INT, MPI_SUM, m_communicator);
#endif
    return diff == 0;
  }

};

template<class scalar, unsigned long N>
inline std::ostream& operator<<(std::ostream& os, OutOfCoreArray<scalar, N> const& obj) { return os << obj.str(); }

}
}  // namespace molpro {

#endif // OUTOFCOREARRAY_H
