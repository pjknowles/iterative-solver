#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_TEMP_FILE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_TEMP_FILE_H
#include <string>
#ifdef HAVE_MPI_H
#include <mpi.h>
#endif
namespace molpro {
namespace linalg {
namespace array {
namespace util {
#ifdef HAVE_MPI_H
/*!
 * @brief Returns random file name for a temporary file
 * @param base_name base name for the file, after which a random character string will follow
 * @param suffix suffix to add to the file name e.g. ".hdf5"
 * @param comm collective over this MPI communicator
 */
std::string temp_file_name(const std::string &base_name, const std::string &suffix, MPI_Comm comm);
#endif
/*!
 * @brief Returns random file name for a temporary file
 * @param base_name base name for the file, after which a random character string will follow
 * @param suffix suffix to add to the file name e.g. ".hdf5"
 */
std::string temp_file_name(const std::string &base_name, const std::string &suffix);

} // namespace util
} // namespace array
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_TEMP_FILE_H
