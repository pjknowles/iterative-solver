#ifndef GCI_SRC_MOLPRO_GCI_ARRAY_ARRAYBASE_H
#define GCI_SRC_MOLPRO_GCI_ARRAY_ARRAYBASE_H

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
 * @warning The distributed array is always open for modification without any locking. We only provide synchronization mechanism
 * to ensure that the RMA calls have completed.
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
  MPI_Comm m_communicator;                            //!< Outer communicator
  std::shared_ptr<molpro::Profiler> m_prof = nullptr; //!< optional profiler
  int m_comm_rank = 0;                                //!< rank in process group
  int m_comm_size = 0;                                //!< size of process group
protected:
  size_t m_dimension = 0; //!< Overall dimension of the array
public:
  ArrayBase() = delete;
  //! Initializes array without allocating any memory
  ArrayBase(size_t dimension, MPI_Comm commun, std::shared_ptr<molpro::Profiler> prof)
      : m_communicator(commun), m_dimension(dimension), m_prof(prof) {}
  ArrayBase(const ArrayBase &source) = default;
  //! Special constructor required by IterativeSolver
  ArrayBase(const ArrayBase &source, int) : ArrayBase(source){};
  [[nodiscard]] ArrayBase &operator=(const ArrayBase &source) = default;
  virtual ~ArrayBase() = default;
  //! Synchronizes all process in this group and ensures any outstanding operations on the array have completed
  virtual void sync() const = 0;
  //! checks if array has been allocated
  virtual bool empty() const = 0;
  //! allocates memory to the array without initializing it with any value. Blocking, collective operation.
  virtual void allocate_buffer() = 0;
  //! Duplicates GA buffer. Requires communicators to be the same. Blocking, collective operation
  virtual void copy_buffer(const ArrayBase &source) = 0;
  //! Provides access to the local portion of the array locking that portion for all other process.
  struct LocalBuffer {
    //! Access local portion of the array
    explicit LocalBuffer(const ArrayBase &source) = 0;
    //! Unlock the buffer for other processes
    virtual ~LocalBuffer() = 0;
    //! Size of the local buffer
    virtual size_t size() const = 0;
    //! Pointer to the start of the buffer
    virtual ArrayBase::value_type *begin() = 0;
    //! Pointer to the end of the buffer (element just after the last one)
    virtual ArrayBase::value_type *end() = 0;
    //! Checks that the current and the other buffers correspond to the same section of their respective arrays
    virtual bool compatible(const LocalBuffer &other) = 0;
    //! Access element at position i relative to begin() without bounds checking
    virtual ArrayBase::value_type &operator[](size_t i) = 0;

    ArrayBase::index_type lo;      //!< first element of local buffer in the array
    ArrayBase::index_type hi;      //!< last element of local buffer in the array (the one before end())
    ArrayBase::value_type *buffer; //!< pointer to the start of the local array buffer
  };
  //! Access the buffer local to this process
  [[nodiscard]] virtual LocalBuffer local_buffer() = 0;
  //! get element at the offset. Blocking.
  [[nodiscard]] virtual ArrayBase::value_type at(size_t offset) const = 0;
  //! Set all local elements to zero.
  virtual void zero(bool with_sync_before = false, bool with_sync_after = false) = 0;
  //! Set all local elements to val.
  virtual void set(ArrayBase::value_type val, bool with_sync_before = false, bool with_sync_after = false) = 0;
  //! Set one element to a scalar. Global operation. @todo rename to put
  virutal void set(size_t ind, ArrayBase::value_type val, bool with_sync_before = false,
                   bool with_sync_after = false) = 0;
  //! Gets buffer[lo:hi] from global array (hi inclusive, i.e. not pythonic). Blocking.
  virtual void get(ArrayBase::index_type lo, ArrayBase::index_type hi,
                   std::vector<ArrayBase::value_type> &buf) const = 0;
  [[nodiscard]] virtual std::vector<ArrayBase::value_type> get(ArrayBase::index_type lo,
                                                               ArrayBase::index_type hi) const = 0;
  //! array[lo:hi] = data[:] (hi inclusive, i.e. not pythonic). Blocking
  virtual void put(ArrayBase::index_type lo, ArrayBase::index_type hi, const ArrayBase::value_type *data) = 0;
  //!  array[lo:hi] += scaling_constant * data[:] (hi inclusive, i.e. not pythonic). Blocking
  virtual void acc(ArrayBase::index_type lo, ArrayBase::index_type hi, const ArrayBase::value_type *data,
                   ArrayBase::value_type scaling_constant = 1.) = 0;
  /*!
   * @brief gets elements with discontinuous indices from array. Blocking
   * @return res[i] = array[indices[i]]
   */
  [[nodiscard]] virtual std::vector<ArrayBase::value_type>
  gather(const std::vector<ArrayBase::index_type> &indices) const = 0;
  /*!
   * @brief array[indices[i]] = data[i]
   * Puts vals of elements with discontinuous indices of array. Blocking.
   */
  [[nodiscard]] virtual void scatter(const std::vector<ArrayBase::index_type> &indices,
                                     const std::vector<ArrayBase::value_type> &data) = 0;
  /*!
   * @brief array[indices[i]] += alpha * vals[i]
   * Accumulates vals of elements into discontinuous indices of array.
   * Atomic, blocking, with on-sided communication
   */
  virtual void scatter_acc(std::vector<ArrayBase::index_type> &indices, const std::vector<ArrayBase::value_type> &data,
                           ArrayBase::value_type alpha) = 0;
  /*!
   * \brief returns n smallest elements
   * Collective operation, must be called by all processes in the group.
   * \param n number of smallest values to be found
   */
  [[nodiscard]] virtual std::list<std::pair<size_t, ArrayBase::value_type>> min_n(size_t n = 1) const = 0;
  /*!
   * \brief returns n largest elements
   * Collective operation, must be called by all processes in the group.
   * \param n number of smallest values to be found
   */
  [[nodiscard]] virtual std::list<std::pair<size_t, ArrayBase::value_type>> max_n(size_t n = 1) const = 0;
  /*!
   * \brief returns n elements that are largest by absolute value
   * Collective operation, must be called by all processes in the group.
   * \param n number of smallest values to be found
   */
  [[nodiscard]] virtual std::list<std::pair<size_t, ArrayBase::value_type>> min_abs_n(size_t n = 1) const = 0;
  /*!
   * \brief returns n elements that are largest by absolute value
   * Collective operation, must be called by all processes in the group.
   * \param n number of smallest values to be found
   */
  [[nodiscard]] virtual std::list<std::pair<size_t, ArrayBase::value_type>> max_abs_n(size_t n = 1) const = 0;
  /*!
   * \brief find the index of n smallest components
   * Collective operation, must be called by all processes in the group.
   * \param n number of smallest values to be found
   * \return offsets in buffer
   */
  [[nodiscard]] virtual std::vector<size_t> minlocN(size_t n = 1) const = 0;
  /*!
   * @brief Copies the whole buffer into a vector. Blocking.
   * @note This is only meant for debugging small arrays!
   */
  [[nodiscard]] virtual std::vector<ArrayBase::value_type> vec() const = 0;
  //! total number of elements, same as overall dimension of array
  virtual size_t size() const = 0;
  //! Checks that arrays are of the same dimensionality
  virtual bool compatible(const ArrayBase &other) const = 0;
  /*!
   * \brief this_array[:] += a * x[:]
   * Add a multiple of another array to this one. Blocking, collective.
   * \param a the pre-factor
   * \param x the other array
   */
  virtual void axpy(ArrayBase::value_type a, const ArrayBase &x, bool with_sync_before = false,
                    bool with_sync_after = false) = 0;
  //! Scale by a constant. Local.
  virtual void scal(ArrayBase::value_type a, bool with_sync_before = false, bool with_sync_after = false) = 0;
  //! Add another array to this. Local
  virtual void add(const ArrayBase &other, bool with_sync_before = false, bool with_sync_after = false) = 0;
  //! Add a constant. Local.
  virtual void add(ArrayBase::value_type a, bool with_sync_before = false, bool with_sync_after = false) = 0;
  //! Subtract another array from this. Local.
  virtual void sub(const ArrayBase &other, bool with_sync_before = false, bool with_sync_after = false) = 0;
  //! Subtract a constant. Local.
  virtual void sub(ArrayBase::value_type a, bool with_sync_before = false, bool with_sync_after = false) = 0;
  //! Take element-wise reciprocal of this. Local. No checks are made for zero values
  virtual void recip(bool with_sync_before = false, bool with_sync_after = false) = 0;
  /*!
   * @brief Scalar product with another array. Collective.
   * Both arrays should be part of the same processor group (same communicator).
   * The result is broadcast to each process.
   */
  virtual ArrayBase::value_type dot(const ArrayBase &other, bool with_sync_before = false) const = 0;
  virtual ArrayBase::value_type dot(const ArrayBase *other, bool with_sync_before = false) const = 0;
  //! this[i] = a[i]*b[i]. Collective.
  virtual void times(const ArrayBase &a, const ArrayBase &b, bool with_sync_before = false,
                     bool with_sync_after = false) = 0;
  /*!
   * \brief this[i] = a[i]/(b[i]+shift). Collective
   * negative? (append? this -=... : this =-...) : (append? this +=... : this =...)
   * \param append Whether to += or =
   * \param negative Whether to scale  right hand side by -1
   */
  virtual void divide(const ArrayBase &a, const ArrayBase *b, ArrayBase::value_type shift = 0, bool append = false,
                      bool negative = false, bool with_sync_before = false, bool with_sync_after = false) = 0;
};

} // namespace molpro::gci::array

#endif // GCI_SRC_MOLPRO_GCI_ARRAY_ARRAYBASE_H
