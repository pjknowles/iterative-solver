#ifndef ITERATIVE_SOLVER_SRC_MOLPRO_LINALG_ARRAY_RESOURCE_H_
#define ITERATIVE_SOLVER_SRC_MOLPRO_LINALG_ARRAY_RESOURCE_H_
#include <filesystem>
#include <molpro/mpi.h>
#include <string>
namespace molpro::linalg::array {

/*!
 * Resources for constructing DistrArray objects
 */
struct Resource {
  MPI_Comm m_mpi_communicator = molpro::mpi::comm_global();
  fs::path m_directory = ".";
  Resource(MPI_Comm communicator = molpro::mpi::comm_global(), fs::path directory = ".")
      : m_mpi_communicator(communicator), m_directory(fs::absolute(fs::path(directory))) {}
};
} // namespace molpro::linalg::array

#endif // ITERATIVE_SOLVER_SRC_MOLPRO_LINALG_ARRAY_RESOURCE_H_
