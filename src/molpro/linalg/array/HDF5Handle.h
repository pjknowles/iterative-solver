#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_HDF5HANDLE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_HDF5HANDLE_H
#include <hdf5.h>
#include <molpro/linalg/array/util/TempHandle.h>
#include <string>
#include <utility>

namespace molpro {
namespace linalg {
namespace array {
namespace util {

/*!
 * @brief Manages opening/closing HDF5 files and groups.
 *
 * Motivation
 * ----------
 * When writing functions or classes that need to write to an hdf5 file there are many ways to pass that information.
 * The user might only have the file name and group name and pass them as strings. Or they could have an hdf5 file
 * object and a group string, or and hdf5 group object, but was the file open with a correct access type or does the
 * group object correspond to a valid file? Does the user want the hdf5 objects to be closed after use?
 *
 * One needs to make multiple constructors and function overloads to ensure the users are kept happy and the procedure
 * might need to be repeated for different classes. The purpose of HDF5Handle is to implement this in one place and let
 * the handle manage opening and closing of hdf5 objects, checking access type and more.
 *
 * Concepts
 * --------
 *   - assignment: handle stores the information about the object
 *   - ownership: handle can open and close the hdf5 objects
 *
 * A file and/or group can be assigned to the handle by passing corresponding file path or group path as a string,
 * or directly providing the relevant hdf5 object.
 * There can only be one file and group assigned to the handle at any one time.
 *
 * Assignment of a file is permanent for each handle, however groups can be re-assigned. Assigning a new
 * group automatically closes the old one, since there can only be one file and group per handle.
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
 * Examples
 * ---------
 * Use a file handle to store arrays in different groups:
 * @code{.cpp}
 * auto h = HDF5Handle(file_name);
 * h.open_file(HDF5Handle::Access::read_write);
 * for (auto & a: arrays){
 *   h.open_group(get_new_group_name()); // reassigning closes the previous group
 *   store_in_new_group(h.group_id(), a);
 * }
 * h.close_file(); // Or let the destructor do it
 * @endcode
 *
 * Use as an argument in a function that needs a group id
 * @code{.cpp}
 * // Without a handle we might need multiple overloads, and perform many checks for validity of the object
 * void write_to_group(std::string file_name, std::string group_name);
 * void write_to_group(hid_t group_id);
 *
 * // Instead we can just take a handle
 * void write_to_group(HDF5Handle& handle){
 *   write_to_group(handle.open_group());
 * }
 * @endcode
 *
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
  explicit HDF5Handle(std::string file);
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
   * @warning If the ownership is not transferred than it is the user's responsibility to close the hdf5 object after
   * the handle is no longer used.
   *
   * @param hid hdf5 id to an open object.
   * @param transfer_ownership whether to take ownership of the object
   */
  explicit HDF5Handle(hid_t hid, bool transfer_ownership = false);

  /*!
   * @brief On destruction closes any hdf5 objects that are owned. If m_erase_on_destroy is true than also erases the
   * underlying file
   */
  virtual ~HDF5Handle();

  /*!
   * @brief Creates a new handle from source. Assignment and ownership of any group or file is copied. Any open and
   * owned objects in source are also opened in the copy.
   * @param source handle to be copied.
   */
  HDF5Handle(const HDF5Handle &source);
  //! Copies source, closing any owned resources and opening a copy of resources owned by source
  HDF5Handle &operator=(const HDF5Handle &source);
  //! Create a handle by transferring ownership from source. Nothing is closed or opened.
  HDF5Handle(HDF5Handle &&source) noexcept;
  //! Close any owned resources and transfer ownership from source
  HDF5Handle &operator=(HDF5Handle &&source) noexcept;

  //! Options for access type when opening files.
  enum class Access { read_only, read_write, none, unknown };
  //!
  /*!
   * @brief Returns access type of the file if it is open or Access::none if it is closed. Can return Access::unknown if
   * there was an error or the access type is not recognised.
   */
  Access access_type() const;
  /*!
   * @brief Open the assigned file.
   *
   * If the file is already open and the access type matches, than it returns the stored hid.
   *
   * The operation fails under the following conditions:
   *   - the named file does not exist and access type is read_only
   *   - the file is already open but with a different access type
   *   - file was assigned without ownership, and was closed on the outside
   *
   * @param type access type
   * @return id of opened hdf5 object or hid_default if operation fails
   */
  virtual hid_t open_file(Access type);
  /*!
   * @brief Open the file passed as input and take ownership of it. Only if a file was not assigned yet.
   *
   * The operation fails under the conditions specified in open_file(Access) and the following:
   *   - a file with a different name was already assigned (reassignment is forbidden, create a new handle)
   *   - a group corresponding to a different file is open (effectively reassignment)
   *
   * @param file name of the file to take ownership of
   * @param type access type
   * @return id of hdf5 object. Returns hid_default if the operation fails or the file was assigned on construction.
   */
  virtual hid_t open_file(const std::string &file, Access type);
  //! Open the group specified on construction, if it is not already open. If the file is not open, opens it with
  //! read_only access
  virtual hid_t open_group();
  /*!
   * @brief Open the group passed as input and take ownership of it.
   *
   * A file must have been already assigned. If the file is not yet open, than it will be opened with read_only access.
   * This will also give the handle ownership of the file object. Since the file is open in read_only state, the
   * handle does not allow modification of the underlying file and there should be no unexpected side-effects.
   *
   * If a group has already been assigned than it will be closed and this handle will be reassigned the new group.
   *
   * @param group name of the group to take ownership of. Should be a full path.
   * @return id of hdf5 object. Returns hid_default if a group was already assigned on construction.
   */
  virtual hid_t open_group(const std::string &group);
  //! Closes the file if it is open, also closes the group. If handle is not owning, the underlying object is not closed
  virtual void close_file();
  //! Closes the group if it is open. If handle is not owning, the underlying object is not closed.
  virtual void close_group();

  //! Returns a copy of hid to the file. The hid is only meaningful if the file is open.
  hid_t file_id() const;
  //! Returns a copy of hid to the group. The hid is only meaningful if the file is open.
  hid_t group_id() const;
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
  bool file_owner() const;
  //! Returns true if the handle owns the group
  bool group_owner() const;
  //! Returns true if the handle does not have an hdf5 object assigned to it
  bool empty() const;

  /*!
   * @brief Sets m_erase_file_on_destroy flag to value.
   * @note Setting the flag to true can fail if the underlying file is not erasable, e.g. handle does not own the file.
   * @param value new value for the flag
   * @return Returns true if the flag was set, or false if it remains the same.
   */
  bool set_erase_file_on_destroy(bool value);
  //! Returns true if the underlying file will be erased when the handle is destroyed
  bool erase_file_on_destroy() const { return m_erase_file_on_destroy; }

  /*!
   * @brief Sets the flag to unlink the group when this object is destroyed.
   * @note The group must own the object.
   * @param value new value for the flag
   * @return Returns true if the flag was set, or false if it remains the same.
   */
  bool set_erase_group_on_destroy(bool value);
  //! Returns true the group will be unlinked on destruction
  bool erase_group_on_destroy() const { return m_erase_group_on_destroy; }

  //! Default value of hid used by the handle. It is always used if the hid is not related to a valid hdf5 object.
  static const hid_t hid_default = -1;

protected:
  hid_t m_file_hid = hid_default;  //!< hdf5 id of the files. Only valid when it is open.
  hid_t m_group_hid = hid_default; //!< hdf5 id of the files. Only valid when it is open.
  std::string m_file_name;         //!< name of the file including full path
  std::string m_group_name;        //!< path to the group in hdf5 file
  bool m_file_owner = false;       //!< flags that the file was open by this instance and should be closed by it
  bool m_group_owner = false;      //!< flags that the group was open by this instance and should be closed by it
  bool m_erase_file_on_destroy =
      false; //!< flags that the underlying file should be erased when the handle is destroyed
  bool m_erase_group_on_destroy = false; //!< flags that the group should be unlinked when the handle is destroyed
  virtual hid_t _open_plist();           //!< returns property list for opening the file
  virtual bool erasable();               //!< returns true if it would be possible to erase the file
};

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

extern template struct TempHandle<HDF5Handle>;
} // namespace util
} // namespace array
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_HDF5HANDLE_H
