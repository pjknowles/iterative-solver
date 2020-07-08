#ifndef GCI_SRC_MOLPRO_GCI_ARRAY_ARRAYBASE_H
#define GCI_SRC_MOLPRO_GCI_ARRAY_ARRAYBASE_H
#include <list>
#include <memory>
#include <mpi.h>
#include <vector>

namespace molpro {
class Profiler;
}

namespace molpro::gci::array {
/*!
 * @brief Declaration of an array distributed across many processes with small set of linear algebra operations
 *
 * This class implements one-sided remote-memory-access operations for getting or putting a copy of any section of the
 * array, provides access to the array data local to the current process, and implements simple linear algebra
 * operations. It also exposes synchronization of processes.
 *
 * This class is designed for the following usage:
 *   - to do simple linear algebra on whole arrays
 *   - getting sections of the array to transform and  accumulate the result into a different array
 *   - initializing the array using put operations
 *   - backing up to and from HDF5 files
 *
 * All RMA operations are blocking, meaning that if putting data into the array the process returns as soon as the data
 * is in the network buffer (not necessarily in the array) and if getting data the process returns when the data
 * is copied into the supplied buffer.
 *
 * The LocalBuffer nested class gives access to the section of distributed array that exists on the current process.
 * There is no locking mechanism and the buffer can still be modified by the RMA operations.
 *
 * The linear algebra operations can be collective or non-collective. In the former, all processes have to participate
 * in the function call, for example dot which requires a collective broadcast. In case of the latter, each process
 * only needs to operate on the local section of the array and no communication is required, for example scaling
 * array by a constant. Collective operations are naturally synchronised.
 *
 * @warning The distributed array is always open for modification without any locking. We only provide synchronization
 * mechanism to ensure that the RMA calls have completed.
 *
 * Example:
 * \code{.cpp}
 *   auto x = Array(n, comm);
 *   // Initialze x
 *   x.allocate();
 *   x.zero(); // non-collective operation. Each process sets local buffer to zero
 *   x.sync(); // Synchronize to make sure all elements are zero
 *   if (rank == 0) { // Simple initialization
 *     x.put(lo, hi, values);
 *     x.scatter(indices, values2);
 *   }
 *   x.sync();
 * \endcode
 *
 */
class ArrayBase {
public:
  using value_type = double;
  using index_type = unsigned long int;
  index_type m_dimension;                             //! number of elements in the array
  MPI_Comm m_communicator;                            //!< Outer communicator
  std::shared_ptr<molpro::Profiler> m_prof = nullptr; //!< optional profiler
  ArrayBase() = delete;
  //! Initializes array without allocating any memory
  ArrayBase(size_t dimension, MPI_Comm commun, std::shared_ptr<molpro::Profiler> prof);
  ArrayBase(const ArrayBase &source) = delete;
  ArrayBase(ArrayBase &&source) = delete;
  ArrayBase &operator=(const ArrayBase &source) = delete;
  ArrayBase &operator=(ArrayBase &&source) = delete;
  virtual ~ArrayBase() = default;

  //! Synchronizes all process in this group and ensures any outstanding operations on the array have completed
  virtual void sync() const;
  //! total number of elements, same as overall dimension of array
  virtual size_t size() const;
  //! Checks that arrays are of the same dimensionality
  virtual bool compatible(const ArrayBase &other) const;
  //! allocates memory to the array without initializing it with any value. Blocking, collective operation.
  virtual void allocate_buffer() = 0;
  //! Duplicates GA buffer. Requires communicators to be the same. Blocking, collective operation
  virtual void copy_buffer(const ArrayBase &source) = 0;
  //! checks if array has been allocated
  virtual bool empty() const;

  /*! @name Local buffer
   * Access the section of the array local to this process
   */
  //! @{
  //! Provides access to the local portion of the array locking that portion for all other process.
  struct LocalBuffer {
    //! Size of the local buffer
    size_t size() { return 1 + hi - lo; };
    //! Pointer to the start of the buffer
    ArrayBase::value_type *begin() { return buffer; };
    //! Pointer to the end of the buffer (element just after the last one)
    ArrayBase::value_type *end() { return buffer + size(); };
    //! Checks that the current and the other buffers correspond to the same section of their respective arrays
    bool compatible(const LocalBuffer &other) { return lo == other.lo && hi == other.hi; };
    //! Access element at position i relative to begin() without bounds checking
    ArrayBase::value_type &at(size_t i) { return buffer[i]; };
    ArrayBase::value_type const &at(size_t i) const { return buffer[i]; };
    const ArrayBase::index_type lo; //!< first element of local buffer in the array
    const ArrayBase::index_type hi; //!< last element of local buffer in the array (the one before end())
  protected:
    ArrayBase::value_type *buffer; //!< pointer to the start of the local array buffer
  };
  //! Access the buffer local to this process
  [[nodiscard]] virtual std::shared_ptr<LocalBuffer> local_buffer() = 0;
  [[nodiscard]] virtual std::shared_ptr<LocalBuffer> local_buffer() const = 0;
  //! @}

  /*! @name One-sided RMA
   * One-sided remote-memory-access operations. They are non-collective
   */
  //! @{
  //! get element at the offset. Blocking.
  [[nodiscard]] virtual value_type at(size_t offset) const = 0;
  //! Set one element to a scalar. Global operation. @todo rename to put
  virtual void set(size_t ind, value_type val) = 0;
  //! Gets buffer[lo:hi] from global array (hi inclusive, i.e. not pythonic). Blocking.
  virtual void get(index_type lo, index_type hi, std::vector<value_type> &buf) const = 0;
  [[nodiscard]] virtual std::vector<value_type> get(index_type lo, index_type hi) const = 0;
  //! array[lo:hi] = data[:] (hi inclusive, i.e. not pythonic). Blocking
  virtual void put(index_type lo, index_type hi, const value_type *data) = 0;
  //!  array[lo:hi] += scaling_constant * data[:] (hi inclusive, i.e. not pythonic). Blocking
  void acc(index_type lo, index_type hi, const value_type *data, value_type scaling_constant = 1.) {
    _acc(lo, hi, data, scaling_constant);
  };
  /*!
   * @brief gets elements with discontinuous indices from array. Blocking
   * @return res[i] = array[indices[i]]
   */
  [[nodiscard]] virtual std::vector<value_type> gather(const std::vector<index_type> &indices) const = 0;
  /*!
   * @brief array[indices[i]] = data[i]
   * Puts vals of elements with discontinuous indices of array. Blocking.
   */
  virtual void scatter(const std::vector<index_type> &indices, const std::vector<value_type> &data) = 0;
  /*!
   * @brief array[indices[i]] += alpha * vals[i]
   * Accumulates vals of elements into discontinuous indices of array.
   * Atomic, blocking, with on-sided communication
   */
  virtual void scatter_acc(std::vector<index_type> &indices, const std::vector<value_type> &data, value_type alpha) = 0;
  /*!
   * @brief Copies the whole buffer into a vector. Blocking.
   * @note This is only meant for debugging small arrays!
   */
  [[nodiscard]] virtual std::vector<value_type> vec() const = 0;
  //! @}

  /*! @name Non-collective
   * Simple linear algebra operations that work on local elements. They are non-collective.
   */
  //! @{
  //! Set all local elements to zero.
  virtual void zero();
  //! Set all local elements to val. @note each process has its own val, there is no communication
  virtual void set(value_type val);
  /*!
   * \brief this_array[:] += a * x[:]
   * Add a multiple of another array to this one. Blocking, collective.
   * \param a the pre-factor
   * \param x the other array
   */
  virtual void axpy(value_type a, const ArrayBase &x);
  //! Scale by a constant. Local.
  virtual void scal(value_type a);
  //! Add another array to this. Local
  virtual void add(const ArrayBase &x);
  //! Add a constant. Local.
  virtual void add(value_type a);
  //! Subtract another array from this. Local.
  virtual void sub(const ArrayBase &x);
  //! Subtract a constant. Local.
  virtual void sub(value_type a);
  //! Take element-wise reciprocal of this. Local. No checks are made for zero values
  virtual void recip();
  //! this[i] *= x[i].
  virtual void times(const ArrayBase &x);
  //! this[i] = x[i]*y[i].
  virtual void times(const ArrayBase &x, const ArrayBase &y);
  //! @}

  /*! @name Collective
   * Collective operations that are blocking and require all processes to call them.
   */
  //!@{
  /*!
   * \brief returns n smallest elements
   * Collective operation, must be called by all processes in the group.
   * \param n number of values to be found
   */
  [[nodiscard]] virtual std::list<std::pair<size_t, value_type>> min_n(size_t n) const;
  /*!
   * \brief returns n largest elements
   * Collective operation, must be called by all processes in the group.
   * \param n number of values to be found
   */
  [[nodiscard]] virtual std::list<std::pair<size_t, value_type>> max_n(size_t n) const;
  /*!
   * \brief returns n elements that are largest by absolute value
   * Collective operation, must be called by all processes in the group.
   * \param n number of smallest values to be found
   */
  [[nodiscard]] virtual std::list<std::pair<size_t, value_type>> min_abs_n(size_t n) const;
  /*!
   * \brief returns n elements that are largest by absolute value
   * Collective operation, must be called by all processes in the group.
   * \param n number of values to be found
   */
  [[nodiscard]] virtual std::list<std::pair<size_t, value_type>> max_abs_n(size_t n) const;
  /*!
   * \brief find the index of n smallest components
   * Collective operation, must be called by all processes in the group.
   * \param n number of values to be found
   * \return offsets in buffer
   */
  [[nodiscard]] virtual std::vector<size_t> min_loc_n(size_t n) const;
  /*!
   * @brief Scalar product with another array. Collective.
   * Both arrays should be part of the same processor group (same communicator).
   * The result is broadcast to each process.
   */
  [[nodiscard]] virtual value_type dot(const ArrayBase &x) const;
  /*!
   * \brief this[i] = x[i]/(y[i]+shift). Collective
   * negative? (append? this -=... : this =-...) : (append? this +=... : this =...)
   * \param x array in the numerator
   * \param y array in the denominator
   * \param append Whether to += or =
   * \param negative Whether to scale  right hand side by -1
   */
  void divide(const ArrayBase &x, const ArrayBase &y, value_type shift = 0, bool append = false,
              bool negative = false) {
    _divide(x, y, shift, append, negative);
  };
  /*! @} */

protected:
  virtual void _acc(index_type lo, index_type hi, const value_type *data, value_type scaling_constant) = 0;
  virtual void _divide(const ArrayBase &x, const ArrayBase &y, value_type shift, bool append, bool negative);
  //! stops application with an error
  virtual void error(const std::string &message) const;
  template <class Compare>[[nodiscard]] std::list<std::pair<index_type, value_type>> extrema(int n) const;
};

} // namespace molpro::gci::array

#endif // GCI_SRC_MOLPRO_GCI_ARRAY_ARRAYBASE_H
