#ifndef STORAGE_H
#define STORAGE_H

#include <cstddef>
#include <fstream>

namespace LinearAlgebra {

/*!
 * \brief The Storage class provides auxiliary storage (usually on an external file) for
 * iterative solvers.
 */
class Storage
{
public:
  /*!
   * \brief Storage
   * \param lengthHint A provided estimate of the total amount of storage that is likely to be needed.
   * \param option An implementation-dependent parameter that controls operation
   */
  Storage(size_t lengthHint=0, int option=0);
  ~Storage();
  /*!
   * \brief Write data to the store.
   * \param buffer Provides the data to be written.
   * \param length Length of data, in bytes.
   * \param address Offset in store, in bytes.
   */
  virtual void write(const char* buffer, size_t length, size_t address);
  /*!
   * \brief Read data from the store.
   * \param buffer Receives the data to be read.
   * \param length Length of data, in bytes.
   * \param address Offset in store, in bytes.
   */
  virtual void read(char *buffer, size_t length, size_t address) const;
  /*!
   * \brief Query the total storage used.
   */
  virtual size_t size() const;
private:
  mutable std::fstream m_file;
  size_t size_; //< total storage (bytes) used
};

}

#endif // STORAGE_H
