#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_HDF5HANDLE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_HDF5HANDLE_H
#include <hdf5.h>
#include <string>
#include <utility>

namespace molpro {
namespace linalg {
namespace array {
namespace util {

/*!
 * @brief Manages opening and closing HDF5 files in serial and parallel modes, opening and closing groups.
 */
class HDF5Handle {
  //! Create a non-owning handle from an already open hdf5 file or group
  HDF5Handle(hid_t hid, bool is_file)
      : m_hid_file{[&hid, &is_file]() -> hid_t { return is_file ? hid : hid_t{}; }()},
        m_hid_group{[&hid, &is_file]() -> hid_t { return !is_file ? hid : hid_t{}; }()} {};

public:
  //! Dummy handle. It should be populated by using open_file() and/or open_group()
  HDF5Handle() = default;
  //! Handle will take ownership of the file, but will not open or close groups.
  HDF5Handle(std::string file) : m_file_name{std::move(file)}, m_file_owner{true} {};
  //! Handle will take ownership of the file and can open/close the specified group
  HDF5Handle(std::string file, std::string group)
      : m_file_name{std::move(file)}, m_group_path{std::move(group)}, m_file_owner{true}, m_group_owner{true} {};
  virtual ~HDF5Handle();

  static HDF5Handle FromFileHandle(hid_t hid) { return HDF5Handle{hid, true}; };
  static HDF5Handle FromGroupHandle(hid_t hid) { return HDF5Handle{hid, false}; };

  //! Open the file specified on construction, if it is not already open
  virtual hid_t open_file();
  //! Open the file passed as input and take ownership of it. Can only be used if no file was specified on construction.
  virtual hid_t open_file(std::string &file);
  //! Open the group specified on construction, if it is not already open. Also opens the file, if it was not open yet.
  virtual hid_t open_group();
  //! Open the group passed as input and take ownership of it. Can only be used if no group was specified on
  //! construction. Also opens the file, if it was not open yet.
  virtual hid_t open_group(std::string &group);
  //! Closes the file if it is open and is owned by this handle
  virtual void close_file();
  //! Closes the group if it is open and is owned by this handle
  virtual void close_group();

  //! Returns a copy of hid to the file. The hid is only meaningful if the file is open.
  hid_t file_id() const { return m_hid_file; }
  //! Returns a copy of hid to the group. The hid is only meaningful if the file is open.
  hid_t group_id() const { return m_hid_group; }
  //! Returns true if the file is open, false otherwise
  bool file_is_open() const;
  //! Returns true if the group is open, false otherwise
  bool group_is_open() const;
  //! Returns true if the handle owns the file
  bool file_owner() const { return m_file_owner; }
  //! Returns true if the handle owns the group
  bool group_owner() const { return m_group_owner; }

protected:
  std::string m_file_name;     //!< name of the file including full path
  std::string m_group_path;    //!< path to the group in hdf5 file
  hid_t m_hid_file = hid_t{};  //!< hdf5 id of the files. Only valid when it is open.
  hid_t m_hid_group = hid_t{}; //!< hdf5 id of the files. Only valid when it is open.
  bool m_file_owner = false;   //!< flags that the file was open by this instance and should be closed by it
  bool m_group_owner = false;  //!< flags that the group was open by this instance and should be closed by it
};

/*!
 * @brief Parallel version of HDF5Handle
 */
class PHDF5Handle : public HDF5Handle {};

//! Returns true if the file exists
bool file_exists(const std::string &fname);

/*!
 * @brief Open an existing hdf5 file in serial mode. File must already exist.
 *
 * @param fname name of the file
 * @param read_only if true open the file with read-only access, otherwise open with read and write access.
 * @return hid of the hdf5 file object. If it is <0, than operation failed.
 */
hid_t hdf5_open_file(const std::string &fname, bool read_only = true);
/*!
 * @brief Open an existing hdf5 file in parallel mode. File must already exist. Collective operation.
 *
 * @param fname name of the file
 * @param communicator mpi communicator. All processes in the communicator must call.
 * @param read_only if true open the file with read-only access, otherwise open with read and write access.
 * @return hid of the hdf5 file object. If it is <0, than operation failed.
 */
hid_t hdf5_open_file(const std::string &fname, MPI_Comm communicator, bool read_only = true);
/*!
 * @brief Open an hdf5 file in serial mode overwriting any existing files under the same name. Opens with read and write
 * access
 *
 * @param fname name of the file
 * @return hid of the hdf5 file object. If it is <0, than operation failed.
 */
hid_t hdf5_create_file(const std::string &fname);
//! Create and open an hdf5 file in parallel mode, overwriting if it already exists. Return its hid, if it is <0, than
//! operation failed.
/*!
 * @brief Open an hdf5 file in parallel mode overwriting any existing files under the same name. Opens with read and
 * write access. Collective operation.
 *
 * @param fname name of the file
 * @param communicator mpi communicator. All processes in the communicator must call.
 * @return hid of the hdf5 file object. If it is <0, than operation failed.
 */
hid_t hdf5_create_file(const std::string &fname, MPI_Comm communicator);

//! Return true if dataset with name dataset_name exists under specified location
bool hdf5_dataset_exists(hid_t location, const std::string &dataset_name);

} // namespace util
} // namespace array
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_HDF5HANDLE_H
