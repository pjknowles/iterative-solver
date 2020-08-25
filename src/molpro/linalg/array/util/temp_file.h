#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_TEMP_FILE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_TEMP_FILE_H
#include <string>
namespace molpro {
namespace linalg {
namespace array {
namespace util {
/*!
 * @brief Returns file name for a temporary file in the current directory that is not already taken
 * @param base_name base name for the file, the suffix will be modified until a valid file name is found
 * @param suffix suffix to add to the file name e.g. ".hdf5"
 */
std::string temp_file_name(const std::string &base_name, const std::string &suffix);

} // namespace util
} // namespace array
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_TEMP_FILE_H
