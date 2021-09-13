#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DISTRARRAYFILE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DISTRARRAYFILE_H

#include "util/fs.h"
#include <molpro/linalg/array/DistrArrayDisk.h>
#include <molpro/mpi.h>

using molpro::mpi::comm_global;

namespace molpro::linalg::array {
/*!
 * @brief Distributed array storing the buffer on disk using temporary local files.
 *
 * On construction, file and corresponding stream object are created and file is then closed.
 * On destruction, file is deleted.
 *
 * @warning Only local operations will be currently supported, if RMA operations are requested, exception will be
 * thrown.
 *
 */
class DistrArrayFile : public DistrArrayDisk {

private:
  fs::path m_directory;
  fs::path m_filename;
  mutable std::unique_ptr<std::fstream>
      m_stream; // use a pointer just in case the implementation of std::fstream is not movable
  mutable std::mutex m_mutex;
  [[nodiscard]] std::tuple<index_type, index_type, index_type>
  local_bounds() const; //! returns local array boundaries and size
public:
  DistrArrayFile() = delete;
  DistrArrayFile(const DistrArrayFile& source);
  //  DistrArrayFile(DistrArrayFile&& source) = default;
  explicit DistrArrayFile(size_t dimension, MPI_Comm comm = comm_global(), const std::string& directory = ".");
  explicit DistrArrayFile(std::unique_ptr<Distribution> distribution, MPI_Comm comm = comm_global(),
                          const std::string& directory = ".");
  explicit DistrArrayFile(const DistrArray& source);

  DistrArrayFile& operator=(const DistrArrayFile& source) = delete;
  DistrArrayFile& operator=(DistrArrayFile&& source) noexcept;

  static DistrArrayFile CreateTempCopy(const DistrArray& source, const std::string& directory = ".");

  friend void swap(DistrArrayFile& x, DistrArrayFile& y) noexcept;

  ~DistrArrayFile() override;

  bool compatible(const DistrArrayFile& source) const;

  //! Dummy
  void erase() override {}
  //! @returns element at given index
  value_type at(index_type ind) const override;
  //! Writes value at a given index
  void set(index_type ind, value_type val) override;
  //! Reads requested values
  void get(index_type lo, index_type hi, value_type* buf) const override;
  //! @returns a vector of values read from file
  std::vector<value_type> get(index_type lo, index_type hi) const override;
  //! Writes requested values to file
  void put(index_type lo, index_type hi, const value_type* data) override;
  //! @warning below functions have to be implemented being pure virtuals, but will return exceptions
  //! taken we don't do collectives here
  void acc(index_type lo, index_type hi, const value_type* data) override;
  std::vector<value_type> gather(const std::vector<index_type>& indices) const override;
  void scatter(const std::vector<index_type>& indices, const std::vector<value_type>& data) override;
  void scatter_acc(std::vector<index_type>& indices, const std::vector<value_type>& data) override;
  std::vector<value_type> vec() const override;
};

} // namespace molpro::linalg::array

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DISTRARRAYFILE_H
