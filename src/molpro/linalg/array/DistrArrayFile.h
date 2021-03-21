#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DISTRARRAYFILE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_DISTRARRAYFILE_H

#include <fstream>
#include <iostream>
#include <set>
 #if ((defined(_MSVC_LANG) && _MSVC_LANG >= 201703L) || (defined(__cplusplus) && __cplusplus >= 201703L)) && defined(__has_include)
 #if __has_include(<filesystem>) && (!defined(__MAC_OS_X_VERSION_MIN_REQUIRED) || __MAC_OS_X_VERSION_MIN_REQUIRED >= 101500) && (!defined(__GNUC__) || __GNUC__ >= 9)
 #define GHC_USE_STD_FS
 #include <filesystem>
 namespace fs = std::filesystem;
 #endif
 #endif
 #ifndef GHC_USE_STD_FS
 #include "ghc/filesystem.h"
 namespace fs = ghc::filesystem;
 #endif

#include "molpro/linalg/array/DistrArrayDisk.h"
#include <molpro/mpi.h>

using molpro::mpi::comm_global;

namespace molpro::linalg::array {
namespace util {
struct FileAttributes {
  explicit FileAttributes(std::fstream stream) : object(std::move(stream)){};
  std::fstream object;
  std::set<std::pair<size_t, size_t>> registry; //! keeps records of start and end lines of an array (object)
};
}
/*!
 * @brief Distributed array storing the buffer on disk using temporary local files.
 *
 * On construction, file and corresponding stream object are being created and file is then closed (to be opened using
 * open_access()). On destruction, file is being deleted.
 *
 * @warning Only local operations will be currently supported, if RMA operations are requested, exception will be
 * thrown.
 *
 */
class DistrArrayFile : public DistrArrayDisk {
private:
  //! creates a file, opens it and @returns m_file fstream
  static std::fstream make_file(const fs::path &dir = fs::current_path());
  std::pair<index_type, index_type> m_frecs; //! start and end indices of an array (object)
  bool m_lrec = false;
  [[nodiscard]] std::tuple<index_type, index_type, index_type> local_bounds() const; //! returns local array boundaries and size
  void update_records(); //! Add a new pair of local array boundaries records; update frecs
public:
  static std::unique_ptr<util::FileAttributes> file;
  //! Constructor for a blank object. The blank is only useful as a temporary. Move a valid object inside the blank to
  //! make it usable.
  DistrArrayFile() = delete;
  DistrArrayFile(const DistrArrayFile &source);
  DistrArrayFile(DistrArrayFile &&source) noexcept;
  explicit DistrArrayFile(size_t dimension, MPI_Comm comm = comm_global(), const std::string &directory = ".");
  explicit DistrArrayFile(std::unique_ptr<Distribution> distribution, MPI_Comm comm = comm_global(), const std::string &directory = ".");
  explicit DistrArrayFile(const DistrArray &source);

  DistrArrayFile &operator=(const DistrArrayFile &source) = delete;
  DistrArrayFile &operator=(DistrArrayFile &&source) noexcept;
  
  static DistrArrayFile CreateTempCopy(const DistrArray &source, const std::string &directory = ".");
  
  friend void swap(DistrArrayFile &x, DistrArrayFile &y) noexcept;

  //! Flushes the buffer if file access is open
  ~DistrArrayFile() override;

  bool compatible(const DistrArrayFile &source) const;

  //! Dummy
  void erase() override {}
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
