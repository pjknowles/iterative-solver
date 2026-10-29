#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DISTRARRAY_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DISTRARRAY_H
#include <cmath>
#include <cstdlib>
#include <list>
#include <map>
#include <memory>
#include <string>
#include <molpro/mpi.h>
#include <vector>
#include <type_traits>
#include <compare>
#include <iterator>
#include <ranges>

#include <molpro/linalg/array/Span.h>
#include <molpro/linalg/array/util/Distribution.h>

namespace molpro::linalg::array {
namespace util {
  /*!
   * Proxy class for APIs that may be used to get or set a single value
   * on the associated array. Depending on whether it is used as a getter
   * or a setter, this wrapper will delegate to the get() and set()
   * functions of the array accordingly.
   * This class is designed to only ever appear as a temporary value that
   * gets immediately converted to the value type or to get destroyed
   * after having written to the value. It is not possible to store
   * instances of this type.
   */
  template<typename Array>
  class ValueProxy {
  public:
    ValueProxy(const ValueProxy &) = delete;
    ValueProxy(ValueProxy &&) = delete;

    // Assigning a proxy writes the *value* it refers to into this proxy's cell,
    // rather than rebinding this proxy to a different array/index -- that keeps
    // *it1 = *it2 (as used by std::copy et al.) working between two DistrArrays.
    ValueProxy &operator=(const ValueProxy &other) { return *this = static_cast<typename Array::value_type>(other); }
    ValueProxy &operator=(ValueProxy &&other) { return *this = static_cast<typename Array::value_type>(other); }
    template <typename OtherArray>
    ValueProxy &operator=(const ValueProxy<OtherArray> &other) {
      return *this = static_cast<typename OtherArray::value_type>(other);
    }

    ValueProxy &operator=(Array::value_type val) {
      m_arr.set(m_idx, std::move(val));

      return *this;
    }

    operator typename Array::value_type() const {
      return m_arr.at(m_idx);
    }

  private:
    std::add_lvalue_reference_t<Array> m_arr;
    Array::index_type m_idx;

    friend std::remove_const_t<Array>;

    ValueProxy(std::add_lvalue_reference_t<Array> arr, Array::index_type idx) : m_arr(arr), m_idx(idx) {}
  };
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
  using index_type = size_t;
  using SparseArray = std::map<index_type, double>;
  using Distribution = util::Distribution<index_type>;

protected:
  index_type m_dimension = 0; //!< number of elements in the array
  MPI_Comm m_communicator;    //!< Outer communicator
  //! Initializes array without allocating any memory
  DistrArray(size_t dimension, MPI_Comm commun);
  DistrArray() = default;

public:
  virtual ~DistrArray() = default;
  // Declared explicitly because the user-declared destructor suppresses the implicit move operations
  DistrArray(const DistrArray&) = default;
  DistrArray(DistrArray&&) noexcept = default;
  DistrArray& operator=(const DistrArray&) = default;
  DistrArray& operator=(DistrArray&&) noexcept = default;

  //! return a copy of the communicator
  MPI_Comm communicator() const { return m_communicator; }
  //! Synchronizes all process in this group and ensures any outstanding operations on the array have completed
  virtual void sync() const;
  //! total number of elements, same as overall dimension of array
  size_t size() const { return m_dimension; };
  //! Checks that arrays are of the same dimensionality
  bool compatible(const DistrArray &other) const;

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
  //! For an array held on disk, which has no local buffer, the number of elements that elementwise operations read
  //! and write at a time; 0 for an array held in memory
  [[nodiscard]] virtual size_t disk_page_size() const { return 0; }
  //! @}

  /*! @name One-sided RMA
   * One-sided remote-memory-access operations. They are non-collective
   */
  //! @{
  //! get element at the offset. Blocking.
  [[nodiscard]] virtual value_type at(index_type ind) const = 0;
  [[nodiscard]] value_type operator[](index_type ind) const;
  [[nodiscard]] util::ValueProxy<DistrArray> operator[](index_type ind);
  //! Set one element to a scalar. Global operation. @todo rename to put
  virtual void set(index_type ind, value_type val) = 0;
  //! Gets buffer[lo:hi) from global array (hi is past-the-end). Blocking.
  virtual void get(index_type lo, index_type hi, value_type *buf) const = 0;
  [[nodiscard]] virtual std::vector<value_type> get(index_type lo, index_type hi) const = 0;
  //! array[lo:hi) = data[:] (hi is past-the-end). Blocking
  virtual void put(index_type lo, index_type hi, const value_type *data) = 0;
  //!  array[lo:hi) += scaling_constant * data[:] (hi is past-the-end). Blocking
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
  //! this[:] = y[:], i.e. the elements of y, and zero elsewhere. Local, each process sets its local section.
  void copy(const SparseArray &y);
  //! @}

  /*! @name Operations with a contiguous array
   * The contiguous array y holds all size() elements on every process, i.e. it is replicated rather than distributed.
   * The distributed array is accessed only in its local section, which is paged through get() and put() if it has no
   * local buffer. y is addressed in place.
   */
  //! @{
  //! this[:] = y[:] over the local section. Local.
  void copy(Span<const value_type> y);
  //! this[:] += a * y[:] over the local section. Local.
  void axpy(value_type a, Span<const value_type> y);
  //! y[:] = this[:]. Collective: each process copies its local section, and then y is replicated over all of them.
  void copy_to(Span<value_type> y) const;
  //! y[:] += a * this[:]. Collective: each process updates its local section, and then y is replicated over all of
  //! them. y must hold the same values on all processes on entry.
  void axpy_to(value_type a, Span<value_type> y) const;
  //! Scalar product <this|y>. Collective, the result is available to each process.
  [[nodiscard]] value_type dot(Span<const value_type> y) const;
  //! The n indices with the largest |this[i]*y[i]|. Collective.
  [[nodiscard]] std::map<size_t, value_type> select_max_dot(size_t n, Span<const value_type> y) const;
  //! @}

  /*! @name Asynchronous linear algebra. No synchronisation on entry or exit.
   */
  //! @{
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

  [[nodiscard]] std::map<size_t, value_type> select_max_dot(size_t n, const DistrArray &y) const;
  [[nodiscard]] std::map<size_t, value_type> select_max_dot(size_t n, const SparseArray &y) const;
  [[nodiscard]] std::map<size_t, value_type> select(size_t n, bool max = false, bool ignore_sign = false) const;
  //! @}

  //! Set all local elements to zero.
  virtual void zero();

  //! stops application with an error
  virtual void error(const std::string &message) const;

protected:
  //! Stops with an error unless the contiguous array y has as many elements as this
  void check_compatible(Span<const value_type> y, const std::string &name) const;
  virtual void _divide(const DistrArray &y, const DistrArray &z, value_type shift, bool append, bool negative);
};

namespace util {
template <typename T, class Compare>
struct CompareAbs {
  constexpr bool operator()(const T &lhs, const T &rhs) const { return Compare()(std::abs(lhs), std::abs(rhs)); }
};
template <class Compare>
[[nodiscard]] std::list<std::pair<DistrArray::index_type, DistrArray::value_type>> extrema(const DistrArray &x, int n);
std::map<size_t, double> select_max_dot_broadcast(size_t n, std::map<size_t, double> &local_selection,
                                                  MPI_Comm communicator);
#ifdef HAVE_MPI_H
//! MPI datatype of DistrArray::index_type; unsigned long is narrower than size_t on Windows
inline MPI_Datatype mpi_index_type() {
  static_assert(sizeof(DistrArray::index_type) == sizeof(unsigned long) ||
                sizeof(DistrArray::index_type) == sizeof(unsigned long long));
  return sizeof(DistrArray::index_type) == sizeof(unsigned long) ? MPI_UNSIGNED_LONG : MPI_UNSIGNED_LONG_LONG;
}
#endif

template<bool is_const>
class DistrArrayIteratorImpl {
  public:
    using value_type = std::conditional_t<
      is_const, DistrArray::value_type, ValueProxy<DistrArray>>;
    using reference = value_type;
    using const_reference = reference;
    using difference_type = std::ptrdiff_t;
    // *this dereferences to a proxy, not an lvalue reference, so there is no meaningful
    // pointer type. iterator_category is spelled out explicitly because std::iterator_traits'
    // legacy (pre-C++20) fallback cannot deduce past input_iterator_tag for a proxy iterator,
    // which would make std::fill/std::copy etc. reject the single-pass-only overloads even
    // though this iterator satisfies std::random_access_iterator (see the static_assert below).
    using pointer = void;
    using iterator_category = std::random_access_iterator_tag;

    DistrArrayIteratorImpl() = default;
    DistrArrayIteratorImpl(std::conditional_t<is_const, const DistrArray, DistrArray> &array,
        DistrArray::index_type pos) : m_arr(&array), m_pos(pos) {}

    DistrArrayIteratorImpl &operator+=(difference_type val) {
      m_pos += val;
      return *this;
    }
    friend DistrArrayIteratorImpl operator+(const DistrArrayIteratorImpl &it, difference_type val) {
      return DistrArrayIteratorImpl(*it.m_arr, it.m_pos + val);
    }
    friend DistrArrayIteratorImpl operator+(difference_type val, const DistrArrayIteratorImpl &it) {
      return it + val;
    }
    DistrArrayIteratorImpl &operator++() {
      ++m_pos;
      return *this;
    }
    DistrArrayIteratorImpl operator++(int) {
      DistrArrayIteratorImpl copy(*this);
      ++m_pos;
      return copy;
    }

    DistrArrayIteratorImpl &operator-=(difference_type val) {
      m_pos -= val;
      return *this;
    }
    friend DistrArrayIteratorImpl operator-(const DistrArrayIteratorImpl &it, difference_type val) {
      return DistrArrayIteratorImpl(*it.m_arr, it.m_pos - val);
    }
    DistrArrayIteratorImpl &operator--() {
      --m_pos;
      return *this;
    }
    DistrArrayIteratorImpl operator--(int) {
      DistrArrayIteratorImpl copy(*this);
      --m_pos;
      return copy;
    }

    reference operator*() const { return (*m_arr)[m_pos]; }

    reference operator[](difference_type offset) const { return (*m_arr)[m_pos + offset]; }

    difference_type operator-(const DistrArrayIteratorImpl &other) const { return m_pos - other.m_pos; }

    bool operator==(const DistrArrayIteratorImpl &) const = default;
    std::strong_ordering operator<=>(const DistrArrayIteratorImpl &) const = default;

  private:
    std::conditional_t<is_const, const DistrArray, DistrArray> *m_arr = nullptr;
    DistrArray::index_type m_pos = 0;
};

static_assert(std::random_access_iterator<DistrArrayIteratorImpl<true>>);
static_assert(std::random_access_iterator<DistrArrayIteratorImpl<false>>);

} // namespace util

using DistrArrayIterator = util::DistrArrayIteratorImpl<false>;
using DistrArrayConstIterator = util::DistrArrayIteratorImpl<true>;

inline DistrArrayIterator begin(DistrArray &array) { return DistrArrayIterator(array, 0); }
inline DistrArrayIterator end(DistrArray &array) { return DistrArrayIterator(array, array.size()); }
inline DistrArrayConstIterator begin(const DistrArray &array) { return DistrArrayConstIterator(array, 0); }
inline DistrArrayConstIterator end(const DistrArray &array) { return DistrArrayConstIterator(array, array.size()); }
inline DistrArrayConstIterator cbegin(const DistrArray &array) { return begin(array); }
inline DistrArrayConstIterator cend(const DistrArray &array) { return end(array); }

} // namespace molpro::linalg::array

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DISTRARRAY_H
