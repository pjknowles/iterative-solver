#ifndef GCI_SRC_MOLPRO_GCI_ARRAY_DISTRARRAYMPI3_H
#define GCI_SRC_MOLPRO_GCI_ARRAY_DISTRARRAYMPI3_H
#include "molpro/gci/array/DistrArray.h"
#include <mpi.h>

namespace molpro {
namespace gci {
namespace array {
namespace util {
class Distribution;
}

/*!
 * @brief Implementation of distributed array using MPI3 RMA operations
 *
 * The array buffer is a window which is open for all processes on creation and
 * remains open until the destruction of the array.
 * @warning Care must be taken that overlapping put and get operations or linear algebra
 * do not cause undefined behaviour.
 */
class DistrArrayMPI3 : public DistrArray {
protected:
protected:
  MPI_Win *m_win = nullptr; //!< window object
  //! distribution of array buffer among processes. Stores start index and size for each
  std::unique_ptr<util::Distribution> m_distribution;
  bool m_allocated = false; //!< whether the window has been created

public:
  DistrArrayMPI3() = delete;
  DistrArrayMPI3(DistrArrayMPI3 &&other) = delete;
  DistrArrayMPI3 &operator=(const DistrArrayMPI3 &&) = delete;

  DistrArrayMPI3(size_t dimension, MPI_Comm commun, std::shared_ptr<Profiler> prof = nullptr);
  //! Copy constructor allocates the buffer if source is not empty
  DistrArrayMPI3(const DistrArray &source);
  DistrArrayMPI3(const DistrArrayMPI3 &source);
  DistrArrayMPI3 &operator=(const DistrArray &source);
  DistrArrayMPI3 &operator=(const DistrArrayMPI3 &source);
  ~DistrArrayMPI3() override;

  void sync() const override;
  void allocate_buffer() override;
  bool empty() const override;

  //! Returns distribution of array buffer among processes with start index and size for each process.
  const util::Distribution &distribution() { return *m_distribution; };

protected:
  struct LocalBufferMPI3 : public DistrArray::LocalBuffer {
    explicit LocalBufferMPI3(DistrArrayMPI3 &source);
  };

public:
  [[nodiscard]] std::shared_ptr<LocalBuffer> local_buffer() override;
  [[nodiscard]] std::shared_ptr<const LocalBuffer> local_buffer() const override;
  [[nodiscard]] value_type at(index_type ind) const override;
  void set(index_type ind, value_type val) override;
  void get(index_type lo, index_type hi, value_type *buf) const override;
  [[nodiscard]] std::vector<value_type> get(index_type lo, index_type hi) const override;
  void put(index_type lo, index_type hi, const value_type *data) override;
  void acc(index_type lo, index_type hi, const value_type *data) override;
  [[nodiscard]] std::vector<value_type> gather(const std::vector<index_type> &indices) const override;
  void scatter(const std::vector<index_type> &indices, const std::vector<value_type> &data) override;
  void scatter_acc(std::vector<index_type> &indices, const std::vector<value_type> &data, value_type alpha) override;
  [[nodiscard]] std::vector<value_type> vec() const override;
  void error(const std::string &message) const override;

protected:
  //! Free the window
  void free_buffer();
  enum class RMAType { get, put, acc, gather, scatter, scatter_acc };
  void _get_put(index_type lo, index_type hi, const value_type *buf, RMAType option);
  //! does gather or scatter or scatter_acc
  void _gather_scatter(const std::vector<index_type> &indices, std::vector<value_type> &data, value_type alpha,
                       RMAType option);
};

} // namespace array
} // namespace gci
} // namespace molpro

#endif // GCI_SRC_MOLPRO_GCI_ARRAY_DISTRARRAYMPI3_H
