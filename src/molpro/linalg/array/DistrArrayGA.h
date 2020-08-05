#ifndef GCI_SRC_MOLPRO_GCI_ARRAY_DISTRARRAYGA_H
#define GCI_SRC_MOLPRO_GCI_ARRAY_DISTRARRAYGA_H
#include <map>

#include "molpro/linalg/array/DistrArray.h"

namespace molpro {
namespace linalg {
namespace array {

/*!
 * @brief Distributed array which uses Global Arrays for managing the array buffer and RMA calls
 */
class DistrArrayGA : public DistrArray {
protected:
  class DistributionGA;
  int m_comm_rank = 0;         //!< rank in process group
  int m_comm_size = 0;         //!< size of process group
  int m_ga_handle = 0;         //!< Global Array handle, needed by GA libary
  int m_ga_pgroup = 0;         //!< Global Array processor group handle
  int m_ga_chunk = 1;          //!< GA chunck size
  bool m_ga_allocated = false; //!< Flags that GA has been allocated
  //! Record every process group created, because GA can only allocate a fixed number of them
  static std::map<MPI_Comm, int> _ga_pgroups;
  std::unique_ptr<DistributionGA> m_distribution = nullptr;

public:
  DistrArrayGA() = default;
  DistrArrayGA(size_t dimension, MPI_Comm commun, std::shared_ptr<Profiler> prof = nullptr);
  //! Make a copy of source, allocating buffer if necessary.
  DistrArrayGA(const DistrArrayGA &source);
  //! Copy source, allocating buffer if necessary.
  DistrArrayGA &operator=(const DistrArrayGA &source);
  //! Take ownership of source leaving it in undefined state. Not collective, no allocation/deallocation
  DistrArrayGA(DistrArrayGA &&source) noexcept;
  //! Take ownership of source leaving it in undefined state. Not collective, no allocation/deallocation
  DistrArrayGA &operator=(DistrArrayGA &&source) noexcept;
  ~DistrArrayGA() override;

  //! swap content of two arrays. Not collective.
  friend void swap(DistrArrayGA &a1, DistrArrayGA &a2);
  void sync() const override;
  void allocate_buffer() override;
  void free_buffer() override;
  bool empty() const override;

protected:
  struct LocalBufferGA : public DistrArray::LocalBuffer {
    explicit LocalBufferGA(DistrArrayGA &source);
  };

  class DistributionGA : public DistrArray::Distribution {
  public:
    DistributionGA() = default;
    DistributionGA(const DistributionGA &) = default;
    DistributionGA(int ga_handle, int n_proc);
    std::pair<int, int> locate_process(index_type lo, index_type hi) const override;
    std::pair<index_type, size_t> range(int process_rank) const override;

  protected:
    int m_ga_handle = 0;
    bool m_dummy = true; //!< marks it as a dummy object, because GA was not allocated yet
    int m_n_proc = 0;    //!< number of processes
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
  //! Checks if index will overflow when converted to int. Calls error() if it does.
  void check_ga_ind_overlow(index_type ind) const;
};

} // namespace array
} // namespace linalg
} // namespace molpro

#endif // GCI_SRC_MOLPRO_GCI_ARRAY_DISTRARRAYGA_H
