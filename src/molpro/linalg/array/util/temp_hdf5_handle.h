#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_TEMP_HDF5_HANDLE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_TEMP_HDF5_HANDLE_H
#include <molpro/linalg/array/HDF5Handle.h>
#include <string>

namespace molpro {
namespace linalg {
namespace array {
namespace util {
/*!
 * @brief Returns a handle to a temporary file that will be erased on its destruction
 *
 * The handle has no group assigned to it.
 *
 * @param base_name base name of the file
 */
HDF5Handle temp_hdf5_handle(const std::string &base_name);

/*!
 * @brief Returns a name for a temporary group that is not already taken
 * @param handle handle with an open file
 * @param base_name base name of the temporary group. It can be absolute path, or relative to group name of handle
 */
std::string temp_group_name(const HDF5Handle &handle, const std::string &base_name);

/*!
 * @brief Returns copy of the handle with a temporary group assigned
 * @param handle handle to modify
 * @param base_name base name of the temporary group. It can be absolute path, or relative to group name of handle
 */
HDF5Handle temp_hdf5_handle_group(const HDF5Handle &handle, const std::string &base_name);

} // namespace util
} // namespace array
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_TEMP_HDF5_HANDLE_H
