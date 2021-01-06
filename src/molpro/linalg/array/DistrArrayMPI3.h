#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DISTRARRAYMPI3_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DISTRARRAYMPI3_H
#include <molpro/linalg/array/DistrArray.h>
#include <mpi.h>

namespace molpro::linalg::array {

/*!
 * @brief Implementation of distributed array using MPI3 RMA operations
 *
 * The array buffer is a window which is open for all processes on creation and
 * remains open until the destruction of the array.
 *
 * @warning Care must be taken that overlapping put and get operations or local buffer modifications
 * do not cause undefined behaviour.
 *
 * @todo I should make implementation more robust by applying properly utilising locks. Current implementation is
 * equivalent to fencing.
 *
 */
class DistrArrayMPI3 : public DistrArray {
protected:
  MPI_Win m_win = MPI_WIN_NULL;                 //!< window object
  std::unique_ptr<Distribution> m_distribution; //!< describes distribution of array among processes
  bool m_allocated = false;                     //!< whether the window has been created

public:
  DistrArrayMPI3() = delete;

  /*!
   * @brief Creates distributed array using MPI3 RMA routines and default distribution.
   *
   * By default array is distributed uniformly with any remainder spread over the first processes so their size is
   * larger by 1.
   *
   * @param dimension overall size of array
   * @param commun MPI communicator
   */
  DistrArrayMPI3(size_t dimension, MPI_Comm commun);

  /*!
   * @brief Creates distributed array using MPI3 RMA routines and specified distribution.
   *
   * @param distribution specification of how array is distributed among processes
   * @param commun MPI communicator
   */
  DistrArrayMPI3(std::unique_ptr<Distribution> distribution, MPI_Comm commun);

  //! Copy constructor allocates the buffer if source is not empty
  DistrArrayMPI3(const DistrArrayMPI3 &source);
  DistrArrayMPI3(const DistrArray &source);
  DistrArrayMPI3(DistrArrayMPI3 &&source) noexcept;
  DistrArrayMPI3 &operator=(const DistrArrayMPI3 &source);
  DistrArrayMPI3 &operator=(DistrArrayMPI3 &&source) noexcept;
  ~DistrArrayMPI3() override;

  friend void swap(DistrArrayMPI3 &a1, DistrArrayMPI3 &a2) noexcept;
  void sync() const override;
  bool empty() const override;

protected:
  struct LocalBufferMPI3 : public DistrArray::LocalBuffer {
    explicit LocalBufferMPI3(DistrArrayMPI3 &source);
  };

public:
  [[nodiscard]] const Distribution &distribution() const override;
  [[nodiscard]] std::unique_ptr<LocalBuffer> local_buffer() override;
  [[nodiscard]] std::unique_ptr<const LocalBuffer> local_buffer() const override;
  [[nodiscard]] value_type at(index_type ind) const override;
  void set(index_type ind, value_type val) override;
  void get(index_type lo, index_type hi, value_type *buf) const override;
  [[nodiscard]] std::vector<value_type> get(index_type lo, index_type hi) const override;
  void put(index_type lo, index_type hi, const value_type *data) override;
  void acc(index_type lo, index_type hi, const value_type *data) override;
  [[nodiscard]] std::vector<value_type> gather(const std::vector<index_type> &indices) const override;
  void scatter(const std::vector<index_type> &indices, const std::vector<value_type> &data) override;
  void scatter_acc(std::vector<index_type> &indices, const std::vector<value_type> &data) override;
  [[nodiscard]] std::vector<value_type> vec() const override;
  void error(const std::string &message) const override;

protected:
  enum class RMAType { get, put, acc, gather, scatter, scatter_acc };
  //! does get or put or accumulate
  void _get_put(index_type lo, index_type hi, const value_type *buf, RMAType option);
  //! does gather or scatter or scatter_acc
  void _gather_scatter(const std::vector<index_type> &indices, std::vector<value_type> &data, RMAType option);
};

} // namespace molpro::linalg::array

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DISTRARRAYMPI3_H
