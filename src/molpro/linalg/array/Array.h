#ifndef GCI_TENSOR_H
#define GCI_TENSOR_H

#include <hdf5.h>
#include <list>
#include <map>
#include <memory>
#include <mpi.h>
#include <string>
#include <vector>

// typedef struct _gci_mpi_comm_dummy *MPI_Comm;
namespace molpro {
class Profiler;
}

namespace molpro::gci::array {
// TODO  1) store the wavefunction on file (hdf5)
//      2) backup current solutions during Davidson
//      3) restart Davidson from previous backup
//      4) use option = 3 in copy constructor to store the wfn on file
//         if set, linear algebra routines should work on local buffers

/*!
 * @brief Wrapper over global array exposing functionality required by IterativeSolver and MixedWavefunction
 */
class Array {
public:
  using value_type = double;
  MPI_Comm m_communicator;                  //!< Outer communicator
  std::shared_ptr<molpro::Profiler> m_prof; //!< optional profiler
  int m_comm_rank;                          //!< rank in process group
  int m_comm_size;                          //!< size of process group
protected:
  size_t m_dimension;  //!< Overall dimension of the direct product Fock space
  int m_ga_handle;     //!< Global Array handle, needed by GA libary
  int m_ga_pgroup;     //!< Global Array processor group handle
  int m_ga_chunk;      //!< GA chunck size
  bool m_ga_allocated; //!< Flags that GA has been allocated
public:
  Array()
      : m_communicator(), m_comm_rank(0), m_comm_size(0), m_dimension(0), m_ga_handle(0), m_ga_pgroup(0), m_ga_chunk(1),
        m_ga_allocated(false) {}

  explicit Array(MPI_Comm comm);

  Array(size_t dimension, MPI_Comm commun, std::shared_ptr<molpro::Profiler> prof = nullptr);

  Array(const Array &source);
  Array(const Array &source, int);

  ~Array();

  bool empty() const; //!< check if GA has been allocated

  void allocate_buffer(); //!< allocates GA buffer. Blocking, collective operation.

  /*!
   * @brief Duplicates GA buffer
   * Requires communicators to be the same. Blocking, collective operation
   */
  void copy_buffer(const Array &source);

  /*!
   * @brief Provides access to the local portion of GA buffer
   */
  class LocalBuffer {
  public:
    explicit LocalBuffer(const Array &source);
    ~LocalBuffer();
    size_t size() const;
    double *begin();
    double *end();
    bool compatible(const LocalBuffer &other);
    double &operator[](size_t i);

    int m_ga_handle, lo, hi, ld;
    double *buffer;
  };

  //! Access the buffer local to this process
  LocalBuffer local_buffer();

  void sync() const;              //!< synchronises all processes in this group and ensures GA operations completed
  double at(size_t offset) const; //!< get element at the offset. Blocking with one-sided communication.
  void zero(bool with_sync_before = false,
            bool with_sync_after = false); ///< Set all local elements to zero
  void set(double val, bool with_sync_before = false,
           bool with_sync_after = false); ///< set all local elements to a scalar
  // Use put instead.
  void set(size_t ind, double val, bool with_sync_before = false,
           bool with_sync_after = false); ///< set one element to a scalar. Global operation

  /*!
   * @brief gets buffer[lo:hi] from global array (hi inclusive, i.e. not pythonic)
   * Blocking with one-sided communication.
   */
  void get(int lo, int hi, std::vector<double> &buf) const;
  std::vector<double> get(int lo, int hi) const;

  /*!
   * @brief  puts data array into GA's buffer[lo:hi] (hi inclusive, i.e. not pythonic)
   * Blocking with one-sided communication
   */
  void put(int lo, int hi, double *data, bool with_fence = false);

  /*!
   * @brief  accumulate data into GA's buffer[lo:hi] (hi inclusive, i.e. not pythonic)
   * Blocking with one-sided communication
   */
  void acc(int lo, int hi, double *buffer, double scaling_constant = 1.);

  /*!
   * @brief gets elements with discontinuos indices from GA
   * Blocking with one-sided communication
   * @return res[i] = ga_buffer[indices[i]]
   */
  std::vector<double> gather(std::vector<int> &indices) const;

  /*!
   * @brief ga_buffer[indices[i]] = vals[i]
   * Puts vals of elements with discontinuos indices into GA.
   * Blocking with one-sided communication.
   */
  void scatter(std::vector<int> &indices, std::vector<double> &vals);

  /*!
   * @brief ga_buffer[indices[i]] += alpha * vals[i]
   * Accumulates vals of elements with discontinuos indices into GA.
   * Atomic, blocking, with on-sided communication
   */
  void scatter_acc(std::vector<int> &indices, std::vector<double> &vals, double alpha);

protected:
  template <class Compare>[[nodiscard]] std::list<std::pair<size_t, double>> extrema(size_t n) const;

public:
  /*!
   * \brief returns n smallest elements
   * Collective operation, must be called by all processes in the group.
   * \param n number of smallest values to be found
   */
  [[nodiscard]] std::list<std::pair<size_t, double>> min_n(size_t n = 1) const;

  /*!
   * \brief returns n largest elements
   * Collective operation, must be called by all processes in the group.
   * \param n number of smallest values to be found
   */
  [[nodiscard]] std::list<std::pair<size_t, double>> max_n(size_t n = 1) const;

  /*!
   * \brief returns n elements that are largest by absolute value
   * Collective operation, must be called by all processes in the group.
   * \param n number of smallest values to be found
   */
  [[nodiscard]] std::list<std::pair<size_t, double>> min_abs_n(size_t n = 1) const;

  /*!
   * \brief returns n elements that are largest by absolute value
   * Collective operation, must be called by all processes in the group.
   * \param n number of smallest values to be found
   */
  [[nodiscard]] std::list<std::pair<size_t, double>> max_abs_n(size_t n = 1) const;

  /*!
   * \brief find the index of n smallest components
   * Collective operation, must be called by all processes in the group.
   * \param n number of smallest values to be found
   * \return offsets in buffer
   */
  std::vector<size_t> minlocN(size_t n = 1) const;

  /*!
   * @brief Copies GA buffer into a local vector
   * Blocking, one-sided communication.
   */
  std::vector<double> vec() const;

  size_t size() const { return m_dimension; } ///< total number of elements

  virtual bool compatible(const Array &other) const; ///< Checks that arrays of the same dimensionality

  /*!
   * \brief axpy Add a multiple of another Wavefunction object to this one
   * Blocking, collective communication
   * \param a the factor defining the multiple
   * \param x the other wavefunction
   * \param with_sync_before synchronise at the start
   * \param with_sync_after synchronise at the end
   */
  void axpy(double a, const Array &x, bool with_sync_before = false, bool with_sync_after = false);
  void axpy(double a, const Array *other, bool with_sync_before = false, bool with_sync_after = false);
  void axpy(double a, const std::map<size_t, double> &x, bool with_sync_before = false, bool with_sync_after = false);

  /*!
   * @brief
   * @param a
   * @param with_sync_before
   * @param with_sync_after
   */
  //! Scale by a constant. Works on local buffer
  void scal(double a, bool with_sync_before = false, bool with_sync_after = false);
  //! Add another array to this. Works on local buffer
  void add(const Array &other, bool with_sync_before = false, bool with_sync_after = false);
  //! Add a constant. Works on local buffer
  void add(double a, bool with_sync_before = false, bool with_sync_after = false);
  //! Subtract another array from this. Works on local buffer
  void sub(const Array &other, bool with_sync_before = false, bool with_sync_after = false);
  //! Subtract a constant. Works on local buffer
  void sub(double a, bool with_sync_before = false, bool with_sync_after = false);
  //! Take element-wise reciprocal of this. Works on local buffer
  void recip(bool with_sync_before = false, bool with_sync_after = false);

  /*!
   * @brief Scalar product with another array.
   * Collective communication. Both arrays should be part of the same processor group (same communicator).
   * The result is broadcast to each process.
   */
  double dot(const Array &other, bool with_sync_before = false) const;
  double dot(const Array *other, bool with_sync_before = false) const;
  double dot(const std::map<size_t, double> &other, bool with_sync_after = false) const;

  Array &operator=(const Array &source) noexcept;
  Array &operator*=(double value);       //!< multiply by a scalar. Collective communication with synchronization
  Array &operator+=(const Array &other); //!< add another array. Collective communication with synchronization
  Array &operator-=(const Array &other); //!< subtract another array. Collective communication with synchronization
  Array &operator+=(double); //!< add a scalar to every element. Collective communication with synchronization
  Array &operator-=(double); //!< subtract a scalar from every element. Collective communication with synchronization
  Array &operator-();        //!< unary minus. Collective communication with synchronization
  //! element-by-element division. Collective communication with synchronization
  Array &operator/=(const Array &other);

  //! this[i] = a[i]*b[i]. Collective communication
  void times(const Array *a, const Array *b, bool with_sync_before = false, bool with_sync_after = false);

  /*!
   * \brief this[i] = a[i]/(b[i]+shift)
   * \param append Whether to do += or =
   * \param negative Whether -shift or +shift
   */
  void divide(const Array *a, const Array *b, double shift = 0, bool append = false, bool negative = false,
              bool with_sync_before = false, bool with_sync_after = false);
}; // class Array

double operator*(const Array &w1, const Array &w2);    ///< inner product of two wavefunctions. Collective communication
Array operator+(const Array &w1, const Array &w2);     ///< add two wavefunctions. Collective communication
Array operator-(const Array &w1, const Array &w2);     ///< subtract two wavefunctions. Collective communication
Array operator/(const Array &w1, const Array &w2);     ///< element-by-element division. Collective communication
Array operator*(const Array &w1, const double &value); ///< multiply by a scalar. Collective communication
Array operator*(const double &value, const Array &w1); ///< multiply by a scalar. Collective communication

} // namespace molpro::gci::array
#endif // GCI_TENSOR_H
