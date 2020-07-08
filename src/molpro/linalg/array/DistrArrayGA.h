#ifndef GCI_SRC_MOLPRO_GCI_ARRAY_DISTRARRAYGA_H
#define GCI_SRC_MOLPRO_GCI_ARRAY_DISTRARRAYGA_H
#include <map>

#include "molpro/gci/array/DistrArray.h"

namespace molpro {
namespace gci {
namespace array {

/*!
 * @brief Distributed array which uses Global Arrays for managing the array buffer and RMA calls
 */
class DistrArrayGA : public DistrArray {
protected:
  int m_comm_rank;     //!< rank in process group
  int m_comm_size;     //!< size of process group
  int m_ga_handle;     //!< Global Array handle, needed by GA libary
  int m_ga_pgroup;     //!< Global Array processor group handle
  int m_ga_chunk;      //!< GA chunck size
  bool m_ga_allocated; //!< Flags that GA has been allocated
  //! Record every process group created, because GA can only allocate a fixed number of them
  static std::map<MPI_Comm, int> _ga_pgroups;

public:
  DistrArrayGA() = delete;
  DistrArrayGA(DistrArrayGA &&other) = delete;
  DistrArrayGA &operator=(const DistrArrayGA &&) = delete;

  DistrArrayGA(size_t dimension, MPI_Comm commun, std::shared_ptr<Profiler> prof = nullptr);
  DistrArrayGA(const DistrArrayGA &other);
  DistrArrayGA &operator=(const DistrArrayGA &);
  ~DistrArrayGA() override;

  void sync() const override;
  void allocate_buffer() override;
  bool empty() const override;

protected:
  struct LocalBufferGA : public DistrArray::LocalBuffer {
    explicit LocalBufferGA(DistrArrayGA &source);
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
  void scatter_acc(std::vector<index_type> &indices, const std::vector<value_type> &data) override;
  [[nodiscard]] std::vector<value_type> vec() const override;
  void error(const std::string &message) const override;

protected:
  //! Checks if index will overflow when converted to int. Calls error() if it does.
  void check_ga_ind_overlow(index_type ind) const;
};

} // namespace array
} // namespace gci
} // namespace molpro

#endif // GCI_SRC_MOLPRO_GCI_ARRAY_DISTRARRAYGA_H
