#ifndef GCI_SRC_MOLPRO_GCI_ARRAY_DISTRARRAY_H
#define GCI_SRC_MOLPRO_GCI_ARRAY_DISTRARRAY_H
#include <cmath>
#include <list>
#include <map>
#include <memory>
#include <mpi.h>
#include <vector>

#include <molpro/linalg/array/Span.h>

namespace molpro {
class Profiler;
}

namespace molpro {
namespace linalg {
namespace array {
namespace util {
template <typename Ind>
class Distribution;
}
/*!
 * @brief Array distributed across many processes supporting remote-memory-access, access to process local buffer, and
 * some linear algebra operations.
 *
 * This class implements one-sided remote-memory-access (RMA) operations for getting or putting a copy of any section of
 * the array, provides access to the array data local to the current process, and implements simple linear algebra
 * operations. It also exposes synchronization of processes, and fencing to ensure RMA operations complete.
 *
 * This class is designed for the following usage:
 *   - to do simple linear algebra on whole arrays
 *   - getting sections of the array to transform and  accumulate the result into a different array
 *   - initializing the array using put operations
 *
 * All RMA operations are blocking, meaning that if putting data into the array the process returns as soon as the data
 * is in the network buffer (not necessarily in the array) and if getting data the process returns when the data
 * is copied into the supplied buffer. Performing synchronisation ensures that all RMA operations complete.
 *
 * The LocalBuffer nested class gives access to the section of distributed array that exists on the current process.
 * It is up to specific implementation of DistributedArray whether exclusive access to the buffer is granted.
 *
 * The linear algebra operations can be collective or non-collective. In the former, all processes have to participate
 * in the function call, for example dot which requires a collective broadcast. In case of the latter, each process
 * only needs to operate on the local section of the array and no communication is required, for example scaling
 * array by a constant. Collective operations are naturally synchronised.
 *
 * @warning The base class does not enforce any locking or exclusive access mechanism. It is up to the specific
 * implementation to decide whether this is necessary. The only rule is that synchronisation call must complete
 * any outstanding RMA and linear algebra operations.
 *
 * Example: blocked matrix vector multiplication
 * \code{.cpp}
 *   // y = A x
 *   auto x = Array(n, comm);
 *   auto y = Array(n, comm);
 *   auto A = Array(n*n, comm);
 *   // Initialze
 *   x.allocate();
 *   x.zero(); // non-collective operation. Each process sets local buffer to zero
 *   y.zero(); // y uses same communicator as x, so does not need a separate synchronisation. It would for RMA.
 *   x.sync(); // Synchronize to make sure all elements are zero
 *   if (rank == 0) { // Simple initialization
 *     x.put(lo, hi, values);
 *     x.scatter(indices, values2);
 *   }
 *   initialize(A); // assume A is stored in row major format
 *   x.sync();
 *   // blocked matrix vector multiplication: y[i] = A[i,j] x[j]
 *   // Let's assume there are nb blocks each of size bs to keep things simple.
 *   // RMA operations need a buffer to copy the data into.
 *   std::vector<double> result_block(bs);
 *   std::vector<double> x_block(bs);
 *   std::vector<double> a_block(bs*bs);
 *   for (auto i = 0; i_col < nb ; ++i){
 *     auto i_lo = i * bs;
 *     auto i_hi = i_lo + bs - 1;
 *     std::fill_n(begin(result_block), bs, 0.);
 *     for (auto j = 0; j < nb ; ++j){
 *       auto j_lo = j * bs;
 *       auto j_hi = j_lo + bs - 1;
 *       if (NextTask()){ // task counter assigns operation to current process
 *         x.get(j_lo, j_hi, x_bloc.data());
 *         A.get((i * nb + j) * bs, (i * nb + j) * bs + bs - 1, a_bloc.data());
 *         // matrix vector multiplication with accumulation into result vector
 *         matrix_vector_multiply(a_bloc, x_bloc, result_block);
 *       }
 *     }
 *     y.accumulate(i_lo, i_hi, result_block.data());
 *   }
 *   y.sync();
 * \endcode
 *
 */
class DistrArray {
public:
  using distributed_array = void; //!< a compile time tag that this is a distributed array
  using value_type = double;
  using index_type = unsigned long int;
  using SparseArray = std::map<unsigned long int, double>;
  using Distribution = util::Distribution<index_type>;

protected:
  index_type m_dimension = 0;              //!< number of elements in the array
  MPI_Comm m_communicator = MPI_COMM_NULL; //!< Outer communicator
  //! Initializes array without allocating any memory
  DistrArray(size_t dimension, MPI_Comm commun, std::shared_ptr<molpro::Profiler> prof);
  DistrArray() = default;

public:
  std::shared_ptr<molpro::Profiler> m_prof = nullptr; //!< optional profiler
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
  //! frees the buffer
  virtual void free_buffer() = 0;
  //! checks if array has been allocated
  virtual bool empty() const;

  /*! @name Local buffer
   * Access the section of the array local to this process
   */
  //! @{
protected:
  //! Provides access to the local portion of the array
  class LocalBuffer : public span::Span<value_type> {
  public:
    virtual ~LocalBuffer() = default;
    using span::Span<value_type>::Span;
    friend void swap(LocalBuffer &, LocalBuffer &) = delete;
    //! Checks that the current and the other buffers correspond to the same section of their respective arrays
    bool compatible(const LocalBuffer &other) const { return start() == other.start() && size() == other.size(); };
    //! Return index to the start of the local buffer section in the distributed array
    size_type start() const { return m_start; }

  protected:
    size_type m_start = 0; //!< index of first element of local buffer in the array
  };

public:
  //! Access the buffer local to this process
  [[nodiscard]] virtual std::unique_ptr<LocalBuffer> local_buffer() = 0;
  [[nodiscard]] virtual std::unique_ptr<const LocalBuffer> local_buffer() const = 0;
  //! Access distribution of the array among processes
  [[nodiscard]] virtual const Distribution &distribution() const = 0;
  //! @}

  /*! @name One-sided RMA
   * One-sided remote-memory-access operations. They are non-collective
   */
  //! @{
  //! get element at the offset. Blocking.
  [[nodiscard]] virtual value_type at(index_type ind) const = 0;
  //! Set one element to a scalar. Global operation. @todo rename to put
  virtual void set(index_type ind, value_type val) = 0;
  //! Gets buffer[lo:hi] from global array (hi inclusive, i.e. not pythonic). Blocking.
  virtual void get(index_type lo, index_type hi, value_type *buf) const = 0;
  [[nodiscard]] virtual std::vector<value_type> get(index_type lo, index_type hi) const = 0;
  //! array[lo:hi] = data[:] (hi inclusive, i.e. not pythonic). Blocking
  virtual void put(index_type lo, index_type hi, const value_type *data) = 0;
  //!  array[lo:hi] += scaling_constant * data[:] (hi inclusive, i.e. not pythonic). Blocking
  virtual void acc(index_type lo, index_type hi, const value_type *data) = 0;
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
   * @brief array[indices[i]] += vals[i]
   * Accumulates vals of elements into discontinuous indices of array.
   * Atomic, blocking, with on-sided communication
   */
  virtual void scatter_acc(std::vector<index_type> &indices, const std::vector<value_type> &data) = 0;
  /*!
   * @brief Copies the whole buffer into a vector. Blocking.
   * @note This is only meant for debugging small arrays!
   */
  [[nodiscard]] virtual std::vector<value_type> vec() const = 0;
  //! @}

  /*! @name Asynchronous linear algebra. No synchronisation on entry or exit.
   */
  //! @{
  //! Set all local elements to val. @note each process has its own val, there is no communication
  virtual void fill(value_type val);
  //! Copies all elements of y. If both arrays are empty than does nothing. If only one is empty, throws an
  //! error.
  virtual void copy(const DistrArray &y);
  /*!
   * @brief Copies elements in a patch of y. If both arrays are empty than does nothing. If only one is empty,
   * throws an error.
   * @param y array to copy
   * @param start index of first element to copy
   * @param end index of last element to copy
   */
  virtual void copy_patch(const DistrArray &y, index_type start, index_type end);
  /*!
   * \brief this[:] += a * y[:]. Throws an error if any array is empty.
   * Add a multiple of another array to this one. Blocking, collective.
   */
  virtual void axpy(value_type a, const DistrArray &y);
  virtual void axpy(value_type a, const SparseArray &y);
  //! Scale by a constant. Local.
  virtual void scal(value_type a);
  //! Add another array to this. Local. Throws error if any array is empty.
  virtual void add(const DistrArray &y);
  //! Add a constant. Local.
  virtual void add(value_type a);
  //! Subtract another array from this. Local. Throws error if any array is empty.
  virtual void sub(const DistrArray &y);
  //! Subtract a constant. Local.
  virtual void sub(value_type a);
  //! Take element-wise reciprocal of this. Local. No checks are made for zero values
  virtual void recip();
  //! this[i] *= y[i]. Throws error if any array is empty.
  virtual void times(const DistrArray &y);
  //! this[i] = y[i]*z[i]. Throws error if any array is empty.
  virtual void times(const DistrArray &y, const DistrArray &z);
  //! @}

  /*! @name Collective linear algebra operations, synchronisation on exit
   */
  //! @{
  /*!
   * @brief Scalar product of two arrays. Collective. Throws error if any array is empty.
   * Both arrays should be part of the same processor group (same communicator).
   * The result is broadcast to each process.
   */
  [[nodiscard]] virtual value_type dot(const DistrArray &y) const;
  [[nodiscard]] virtual value_type dot(const SparseArray &y) const;

  /*!
   * @brief this[i] = y[i]/(z[i]+shift). Collective. Throws error if any array is empty.
   * @code{.cpp}
   * negative? (append? this -=... : this =-...) : (append? this +=... : this =...)
   * @endcode
   * @param y array in the numerator
   * @param z array in the denominator
   * @param shift denominator shift
   * @param append Whether to += or =
   * @param negative Whether to scale  right hand side by -1
   */
  void divide(const DistrArray &y, const DistrArray &z, value_type shift = 0, bool append = false,
              bool negative = false) {
    _divide(y, z, shift, append, negative);
  }

  /*!
   * @brief returns n smallest elements in array x
   * Collective operation, must be called by all processes in the group.
   * @return list of index and value pairs, or empty list if array is empty.
   */
  [[nodiscard]] std::list<std::pair<index_type, value_type>> min_n(int n) const;

  /*!
   * \brief returns n largest elements in array x
   * Collective operation, must be called by all processes in the group.
   * @return list of index and value pairs, or empty list if array is empty.
   */
  [[nodiscard]] std::list<std::pair<index_type, value_type>> max_n(int n) const;

  /*!
   * \brief returns n elements that are largest by absolute value in array x
   * Collective operation, must be called by all processes in the group.
   * @return list of index and value pairs, or empty list if array is empty.
   */
  [[nodiscard]] std::list<std::pair<index_type, value_type>> min_abs_n(int n) const;

  /*!
   * \brief returns n elements that are largest by absolute value in array x
   * Collective operation, must be called by all processes in the group.
   * @return list of index and value pairs, or empty list if array is empty.
   */
  [[nodiscard]] std::list<std::pair<index_type, value_type>> max_abs_n(int n) const;

  /*!
   * \brief find the index of n smallest components in array x
   * Collective operation, must be called by all processes in the group.
   * @return list of indices for smallest n values, or empty list if array is empty.
   */
  [[nodiscard]] std::vector<index_type> min_loc_n(int n) const;
  //! @}

  //! Set all local elements to zero.
  virtual void zero();

  //! stops application with an error
  virtual void error(const std::string &message) const;

protected:
  virtual void _divide(const DistrArray &y, const DistrArray &z, value_type shift, bool append, bool negative);
};

namespace util {
template <typename T, class Compare>
struct CompareAbs {
  constexpr bool operator()(const T &lhs, const T &rhs) const { return Compare()(std::abs(lhs), std::abs(rhs)); }
};
template <class Compare>
[[nodiscard]] std::list<std::pair<DistrArray::index_type, DistrArray::value_type>> extrema(const DistrArray &x, int n);
} // namespace util
} // namespace array
} // namespace linalg
} // namespace molpro

#endif // GCI_SRC_MOLPRO_GCI_ARRAY_DISTRARRAY_H
