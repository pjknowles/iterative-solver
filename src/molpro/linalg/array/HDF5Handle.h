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
 * @brief Manages opening/closing HDF5 files and groups.
 *
 * A file and/or group can be assigned to the handle by passing corresponding file path or group path as a string,
 * or directly providing the relevant hdf5 object.
 *
 * If the assignment is done by passing strings than the handle takes ownership of the hdf5 objects and
 * can open and close the file and/or group.
 * Opening is done by calling open_file()/open_group(), and closing can be done either through
 * close_file()/close_group() or by the destructor.
 *
 * The assignment can be done by passing the hid in which case the user can transfer ownership of the underlying
 * hdf5 object allowing handle to close it, or they can keep ownership to themselves thus keeping responsibility
 * for closing the object when it is not needed.
 *
 * Concepts
 * --------
 *   - assignment: handle stores the information about the object
 *   - ownership: handle can open and close the hdf5 objects
 *
 */
class HDF5Handle {
public:
  /*!
   * @brief Create a dummy handle taking on default values, with no underlying assigned objects.
   *
   * It should be populated by using open_file() and/or open_group()
   */
  HDF5Handle() = default;
  //! Create a handle with file assigned to it with ownership
  HDF5Handle(std::string file);
  //! Create a handle with file and group assigned to it with ownership
  HDF5Handle(std::string file, std::string group);
  /*!
   * @brief Create a handle from an already open hdf5 file or group with option to transfer ownership
   *
   * Only one hdf5 object can be assigned. If it is a file object, than later calls to open_group() can assign
   * a group with ownership making sure that it can be closed.
   * Assigning a group object implies assignment of the overlying file with the same ownership.
   * If the ownership is not transferred than the overlying hdf5 file object must remain open.
   * Furthermore, opening a file that does not contain a group object is erroneous.
   *
   * File name containing the hdf5 object is always stored, and if it is a group object than group name is also stored.
   *
   * The handle will be in a dummy state in the following cases:
   *   - hid object is not valid (see H5Iis_valid)
   *   - hid object does not correspond to a file or a group
   *
   * @warning If the ownership is not transferred than it is the user's responsibility to close the file after the
   * handle is no longer used.
   *
   * @param hid hdf5 id to an open object.
   * @param transfer_ownership whether to take ownership of the object
   */
  HDF5Handle(hid_t hid, bool transfer_ownership = false);

  virtual ~HDF5Handle();

  //! Creates a new handle from source. If source is an owner of file or group, than they are opened again.
  HDF5Handle(const HDF5Handle &source);
  //! Copies source, closing any owned resources and opening a copy of resources owned by source
  HDF5Handle &operator=(const HDF5Handle &source);
  //! Create a handle by transferring ownership from source
  HDF5Handle(HDF5Handle &&source);
  //! Close any owned resources and transfer ownership from source
  HDF5Handle &operator=(HDF5Handle &&source);

  //! Options for access type when opening files.
  enum class Access { read_only, read_write };
  /*!
   * @brief Open the assigned file.
   *
   * If the file is already open and the access type matches, than it returns the stored hid.
   *
   * The operation fails under the following conditions:
   *   - the named file does not exist and access type is read_only
   *   - the file is already open but with a different access type
   *   - file was assigned without ownership TODO unit test
   *
   * @param type access type
   * @return id of opened hdf5 object or hid_default if operation fails
   */
  virtual hid_t open_file(Access type);
  /*!
   * @brief Open the file passed as input and take ownership of it. Only if a file was not assigned yet.
   *
   * The operation fails under the following conditions:
   *   - the file is already open but with a different access type
   *   - a file with a different name was already assigned
   *   - the named file does not exist and access type is read_only
   *   - a group corresponding to a different file is open TODO unit test
   *
   * @param file name of the file to take ownership of
   * @return id of hdf5 object. Returns hid_default if the operation fails or the file was assigned on construction.
   */
  virtual hid_t open_file(const std::string &file, Access type);
  //! Open the group specified on construction, if it is not already open. If the file is not open, opens it with
  //! read_only access
  virtual hid_t open_group();
  /*!
   * @brief Open the group passed as input and take ownership of it.
   *
   * A file must have been assigned already, but not a group. Opens the file with read_only access, if it was not open
   * yet.
   *
   * @param group name of the group to take ownership of. Should be a full
   * @return id of hdf5 object. Returns hid_default if a group was already assigned on construction.
   */
  virtual hid_t open_group(std::string &group);
  //! Closes the file if it is open and is owned by this handle
  virtual void close_file();
  //! Closes the group if it is open and is owned by this handle
  virtual void close_group();

  //! Returns a copy of hid to the file. The hid is only meaningful if the file is open.
  hid_t file_id() const { return m_file_hid; }
  //! Returns a copy of hid to the group. The hid is only meaningful if the file is open.
  hid_t group_id() const { return m_group_hid; }
  //! Returns the file name of an object assigned to the handle, even if it is not owned by it.
  std::string file_name() const;
  //! Returns the group name of an object assigned to the handle, even if it is not owned by it.
  std::string group_name() const;
  /*!
   * @brief Checks that the file handle is open, even if it is not owned.
   *
   * File is considered open if the underlying hdf5 object is valid (see H5Iis_valid)
   *
   * @return true if the file is open and can be used with hdf5 operations, false otherwise
   */
  bool file_is_open() const;
  /*!
   * @brief Checks that the group handle is open, even if it is not owned.
   *
   * Group is considered open if the underlying hdf5 object is valid (see H5Iis_valid)
   *
   * @return true if the group is open and can be used with hdf5 operations, false otherwise
   */
  bool group_is_open() const;
  //! Returns true if the handle owns the file
  bool file_owner() const { return m_file_owner; }
  //! Returns true if the handle owns the group
  bool group_owner() const { return m_group_owner; }
  //! Returns true if the handle does not have an hdf5 object assigned to it
  bool empty() const { return m_file_hid == hid_default && m_group_hid == hid_default; }

  //! Default value of hid used by the handle. It is always used if the hid is not related to a valid hdf5 object.
  static const hid_t hid_default = -1;

protected:
  hid_t m_file_hid = hid_default;  //!< hdf5 id of the files. Only valid when it is open.
  hid_t m_group_hid = hid_default; //!< hdf5 id of the files. Only valid when it is open.
  std::string m_file_name;         //!< name of the file including full path
  std::string m_group_name;        //!< path to the group in hdf5 file
  bool m_file_owner = false;       //!< flags that the file was open by this instance and should be closed by it
  bool m_group_owner = false;      //!< flags that the group was open by this instance and should be closed by it
  bool m_file_is_open = false;     //!< flags that the file was opened by the class
  bool m_group_is_open = false;    //!< flags that the group was opened by the class
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
//! Returns the name of hdf5 object. See H5Iget_name for valid object ids.
std::string hdf5_get_object_name(hid_t id);
//! Returns file name of hdf5 object. If the input object is not a file than returns name of while to which it belongs.
//! See H5Fget_name
std::string hdf5_get_file_name(hid_t id);
/*!
 * @brief Starting from id object checks that each link in the path exists.
 *
 * E.g. if id is file and path="/g1/g2/g3" it checks that "/g1", "/g1/g2", and "/g1/g2/g3" all exist.
 *
 * @param id hdf5 object id
 * @param path path to the link relative to the id. If id is a file object, should include root "/"
 * @return >0 if link exists, 0 if it doesnt, <0 if it failed or path is empty.
 */
//!
htri_t hdf5_link_exists(hid_t id, std::string path);

} // namespace util
} // namespace array
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_HDF5HANDLE_H
