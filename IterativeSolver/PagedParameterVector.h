#ifndef PAGEDPARAMETERVECTOR_H
#define PAGEDPARAMETERVECTOR_H
#include "ParameterVector.h"
#include "Storage.h"
#include <fstream>

namespace IterativeSolver {

  /*!
   * \brief A class that implements ParameterVector with data held on backing store
   */
  class PagedParameterVector : public ParameterVector
  {
  public:
    /*!
   * \brief Construct an object without any data.
   */
    PagedParameterVector(size_t length=0);
    PagedParameterVector(const PagedParameterVector& source);
    ~PagedParameterVector();
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
    PagedParameterVector& operator=(const PagedParameterVector& other);

    // Every child of ParameterVector needs exactly this
    PagedParameterVector* clone() const { return new PagedParameterVector(*this); }

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
    std::vector<ParameterScalar> m_cache;
    bool m_cacheDirty;
    long long m_cacheOffset;
    void write(ParameterScalar* const buffer, size_t length, size_t offset);
    void read(ParameterScalar* buffer, size_t length, size_t offset) const;
  public:
    void put(ParameterScalar* const buffer, size_t length, size_t offset);
    void get(ParameterScalar* buffer, size_t length, size_t offset) const;
    size_t size() const;
    std::string str() const;
  };

}

#endif // PAGEDPARAMETERVECTOR_H
