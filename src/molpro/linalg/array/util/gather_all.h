#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_GATHER_ALL_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_GATHER_ALL_H
#include <molpro/linalg/array/util/Distribution.h>
#include <mpi.h>

namespace molpro::linalg::array::util {

/*!
 * @brief Replicate data of a full container on all the processes based on distributed pieces
 *
 * @param distr distribution of the container, which local data chuncks are to be replicated
 * @param commun MPI communicator
 * @param first_elem pointer to the beginning of the container, which data to be replicated
 */
void gather_all(const Distribution<size_t>& distr, MPI_Comm commun, double* first_elem) {
  int nproc, mpi_rank;
  MPI_Comm_size(commun, &nproc);
  MPI_Comm_rank(commun, &mpi_rank);
  int chunks[nproc], displs[nproc];
  for (int i = 0; i < nproc; i++) {
    displs[i] = distr.range(i).first;
    chunks[i] = distr.range(i).second - distr.range(i).first;
  }
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, first_elem, chunks, displs, MPI_DOUBLE, commun);
}

} // namespace molpro::linalg::array::util

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_GATHER_ALL_H
