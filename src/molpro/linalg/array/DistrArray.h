#ifndef GCI_SRC_MOLPRO_GCI_ARRAY_DISTRARRAY_H
#define GCI_SRC_MOLPRO_GCI_ARRAY_DISTRARRAY_H
#include <list>
#include <memory>
#include <mpi.h>
#include <vector>

namespace molpro {
class Profiler;
}

namespace molpro::gci::array {
/*!
 * @brief Array distributed across many processes with small set of linear algebra operations
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
class DistrArray {
public:
  using value_type = double;
  using index_type = unsigned long int;

protected:
  index_type m_dimension;  //! number of elements in the array
  MPI_Comm m_communicator; //!< Outer communicator
public:
  std::shared_ptr<molpro::Profiler> m_prof = nullptr; //!< optional profiler
  DistrArray() = delete;
  //! Initializes array without allocating any memory
  DistrArray(size_t dimension, MPI_Comm commun, std::shared_ptr<molpro::Profiler> prof);
  DistrArray(const DistrArray &source) = delete;
  DistrArray(DistrArray &&source) = delete;
  DistrArray &operator=(const DistrArray &source) = delete;
  DistrArray &operator=(DistrArray &&source) = delete;
  virtual ~DistrArray() = default;

  //! return a copy of the communicator
  MPI_Comm communicator() const { return m_communicator; }
  //! Synchronizes all process in this group and ensures any outstanding operations on the array have completed
  virtual void sync() const;
  //! total number of elements, same as overall dimension of array
  size_t size() const { return m_dimension; };
  //! Checks that arrays are of the same dimensionality
  bool compatible(const DistrArray &other) const;
  //! allocates memory to the array without initializing it with any value. Blocking, collective operation.
  virtual void allocate_buffer() = 0;
  //! checks if array has been allocated
  virtual bool empty() const;

  /*! @name Local buffer
   * Access the section of the array local to this process
   */
  //! @{
protected:
  //! Provides access to the local portion of the array locking that portion for all other process.
  class LocalBuffer {
  public:
    //! Size of the local buffer
    size_t size() const { return hi - lo; };
    //! Pointer to the start of the buffer
    DistrArray::value_type *begin() { return buffer; };
    const DistrArray::value_type *begin() const { return buffer; };
    const DistrArray::value_type *cbegin() const { return begin(); };
    //! Pointer to the end of the buffer (element just after the last one)
    DistrArray::value_type *end() { return buffer + size(); };
    const DistrArray::value_type *end() const { return buffer + size(); };
    const DistrArray::value_type *cend() const { return end(); };
    //! Checks that the current and the other buffers correspond to the same section of their respective arrays
    bool compatible(const LocalBuffer &other) const { return lo == other.lo && hi == other.hi; };
    //! Access element at position i relative to begin() without bounds checking
    DistrArray::value_type &at(size_t i) { return buffer[i]; };
    DistrArray::value_type const &at(size_t i) const { return buffer[i]; };
    DistrArray::index_type lo;      //!< index of first element of local buffer in the array
    DistrArray::index_type hi;      //!< index of one past the last element in the buffer (same as end())
    DistrArray::value_type *buffer; //!< pointer to the start of the local array buffer
  };

public:
  //! Access the buffer local to this process
  [[nodiscard]] virtual std::shared_ptr<LocalBuffer> local_buffer() = 0;
  [[nodiscard]] virtual std::shared_ptr<const LocalBuffer> local_buffer() const = 0;
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

  //! Set all local elements to zero.
  virtual void zero();

  //! stops application with an error
  virtual void error(const std::string &message) const;

protected:
  virtual void _acc(index_type lo, index_type hi, const value_type *data, value_type scaling_constant) = 0;
};

//! Set all local elements of array x to val. @note each process has its own val, there is no communication
void fill(DistrArray &x, DistrArray::value_type val);
//! Copies all elements of y into x. If both arrays are empty than does nothing. If only one is empty, throws an error.
void copy(DistrArray &x, const DistrArray &y);
/*!
 * @brief Copies elements in a patch of y into x. If both arrays are empty than does nothing. If only one is empty,
 * throws an error.
 * @param x array to be copied into
 * @param y array to copy
 * @param start index of first element to copy
 * @param end index of last element to copy
 */
void copy_patch(DistrArray &x, const DistrArray &y, DistrArray::index_type start, DistrArray::index_type end);
/*!
 * \brief x[:] += a * y[:]. Throws error if any array is empty.
 * Add a multiple of another array to this one. Blocking, collective.
 */
void axpy(DistrArray &x, DistrArray::value_type a, const DistrArray &y);
//! Scale by a constant. Local.
void scal(DistrArray &x, DistrArray::value_type a);
//! Add another array to this. Local. Throws error if any array is empty.
void add(DistrArray &x, const DistrArray &y);
//! Add a constant. Local.
void add(DistrArray &x, DistrArray::value_type a);
//! Subtract another array from this. Local. Throws error if any array is empty.
void sub(DistrArray &x, const DistrArray &y);
//! Subtract a constant. Local.
void sub(DistrArray &x, DistrArray::value_type a);
//! Take element-wise reciprocal of this. Local. No checks are made for zero values
void recip(DistrArray &x);
//! x[i] *= y[i]. Throws error if any array is empty.
void times(DistrArray &x, const DistrArray &y);
//! x[i] = y[i]*z[i]. Throws error if any array is empty.
void times(DistrArray &x, const DistrArray &y, const DistrArray &z);

/*!
 * @brief Scalar product of two arrays. Collective. Throws error if any array is empty.
 * Both arrays should be part of the same processor group (same communicator).
 * The result is broadcast to each process.
 */
[[nodiscard]] DistrArray::value_type dot(const DistrArray &x, const DistrArray &y);

/*!
 * @brief x[i] = y[i]/(z[i]+shift). Collective. Throws error if any array is empty.
 * @code{.cpp}
 * negative? (append? this -=... : this =-...) : (append? this +=... : this =...)
 * @endcode
 * @param x result array
 * @param y array in the numerator
 * @param z array in the denominator
 * @param shift denominator shift
 * @param append Whether to += or =
 * @param negative Whether to scale  right hand side by -1
 */
void divide(DistrArray &x, const DistrArray &y, const DistrArray &z, DistrArray::value_type shift = 0,
            bool append = false, bool negative = false);

/*!
 * @brief returns n smallest elements in array x
 * Collective operation, must be called by all processes in the group.
 * @return list of index and value pairs, or empty list if array is empty.
 */
[[nodiscard]] std::list<std::pair<DistrArray::index_type, DistrArray::value_type>> min_n(const DistrArray &x, int n);

/*!
 * \brief returns n largest elements in array x
 * Collective operation, must be called by all processes in the group.
 * @return list of index and value pairs, or empty list if array is empty.
 */
[[nodiscard]] std::list<std::pair<DistrArray::index_type, DistrArray::value_type>> max_n(const DistrArray &x, int n);

/*!
 * \brief returns n elements that are largest by absolute value in array x
 * Collective operation, must be called by all processes in the group.
 * @return list of index and value pairs, or empty list if array is empty.
 */
[[nodiscard]] std::list<std::pair<DistrArray::index_type, DistrArray::value_type>> min_abs_n(const DistrArray &x,
                                                                                             int n);

/*!
 * \brief returns n elements that are largest by absolute value in array x
 * Collective operation, must be called by all processes in the group.
 * @return list of index and value pairs, or empty list if array is empty.
 */
[[nodiscard]] std::list<std::pair<DistrArray::index_type, DistrArray::value_type>> max_abs_n(const DistrArray &x,
                                                                                             int n);

/*!
 * \brief find the index of n smallest components in array x
 * Collective operation, must be called by all processes in the group.
 * @return list of indices for smallest n values, or empty list if array is empty.
 */
[[nodiscard]] std::vector<DistrArray::index_type> min_loc_n(const DistrArray &x, int n);

namespace util {
template <class Compare>
[[nodiscard]] std::list<std::pair<DistrArray::index_type, DistrArray::value_type>> extrema(const DistrArray &x, int n);
}

} // namespace molpro::gci::array

#endif // GCI_SRC_MOLPRO_GCI_ARRAY_DISTRARRAY_H
