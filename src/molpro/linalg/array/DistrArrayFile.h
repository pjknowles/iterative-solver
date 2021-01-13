#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DISTRARRAYFILE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DISTRARRAYFILE_H

#include <fstream>
#include <iostream>
#ifdef USE_EXPERIMENTAL_FILESYSTEM
#include <experimental/filesystem>
#else
#include <filesystem>
#endif

#include "molpro/linalg/array/DistrArrayDisk.h"

namespace molpro::linalg::array {

#ifdef USE_EXPERIMENTAL_FILESYSTEM
namespace fs = std::experimental::filesystem;
#else
namespace fs = std::filesystem;
#endif
/*!
 * @brief Distributed array storing the buffer on disk using temporary local files.
 *
 * On construction, file and corresponding stream object are being created and file is then closed (to be opened using
 * open_access()). On destruction, file is being deleted.
 *
 * @warning Only local operations will be currently supported, if RMA operations are requested, exception will be thrown.
 *
 */
class DistrArrayFile : public DistrArrayDisk {
protected:
  fs::path m_dir = fs::current_path();
  mutable std::fstream m_file;
  //! creates a file, opens it and @returns m_file fstream
  std::fstream make_file();
public:
  //! Constructor for a blank object. The blank is only useful as a temporary. Move a valid object inside the blank to
  //! make it usable.
  DistrArrayFile();
  DistrArrayFile(const DistrArrayFile &source);
  DistrArrayFile(DistrArrayFile &&source) noexcept;
  explicit DistrArrayFile(size_t dimension, MPI_Comm comm = MPI_COMM_WORLD, const std::string &directory = ".");
  explicit DistrArrayFile(std::unique_ptr<Distribution> distribution, MPI_Comm comm = MPI_COMM_WORLD, const std::string &directory = ".");
  explicit DistrArrayFile(const DistrArray &source);
  
  DistrArrayFile &operator=(const DistrArrayFile &source) = delete;
  DistrArrayFile &operator=(DistrArrayFile &&source) noexcept;
  
  static DistrArrayFile CreateTempCopy(const DistrArray &source, const std::string &directory = ".");
  
  friend void swap(DistrArrayFile &x, DistrArrayFile &y) noexcept;
  
  //! Flushes the buffer if file access is open
  ~DistrArrayFile() override;
  
  bool compatible(const DistrArrayFile &source) const;
 
  //! Dummy
  void open_access() override;
  //! Dummy
  void close_access() override;
  //! @returns true if array is not accessible through file nor memory view. Returns false otherwise.
  bool empty() const override;
  //! Dummy
  void erase() override;
  //! @returns element at given index
  value_type at(index_type ind) const override;
  //! Writes value at a given index
  void set(index_type ind, value_type val) override;
  //! Reads requested values
  void get(index_type lo, index_type hi, value_type *buf) const override;
  //! @returns a vector of values read from file
  std::vector<value_type> get(index_type lo, index_type hi) const override;
  //! Writes requested values to file
  void put(index_type lo, index_type hi, const value_type *data) override;
  //! @warning below functions have to be implemented being pure virtuals, but will return exceptions
  //! taken we don't do collectives here
  void acc(index_type lo, index_type hi, const value_type *data) override;
  std::vector<value_type> gather(const std::vector<index_type> &indices) const override;
  void scatter(const std::vector<index_type> &indices, const std::vector<value_type> &data) override;
  void scatter_acc(std::vector<index_type> &indices, const std::vector<value_type> &data) override;
  std::vector<value_type> vec() const override;
};

} // namespace molpro::linalg::array

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DISTRARRAYFILE_H
