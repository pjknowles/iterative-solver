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
HDF5Handle temp_hdf5_handle(std::string base_name);

} // namespace util
} // namespace array
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_TEMP_HDF5_HANDLE_H
