#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_GATHER_ALL_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_GATHER_ALL_H
#include <molpro/linalg/array/util/Distribution.h>
#include <molpro/mpi.h>
#include <vector>

namespace molpro::linalg::array::util {

/*!
 * @brief Replicate data of a full container on all the processes based on distributed pieces
 *
 * @param distr distribution of the container, which local data chuncks are to be replicated
 * @param commun MPI communicator
 * @param first_elem pointer to the beginning of the container, which data to be replicated
 *
 * Does nothing without MPI, or on a single process, when the container already holds all the data.
 */
inline void gather_all(const Distribution<size_t>& distr, MPI_Comm commun, double* first_elem) {
#ifdef HAVE_MPI_H
  int nproc;
  MPI_Comm_size(commun, &nproc);
  if (nproc == 1)
    return;
  std::vector<int> chunks(nproc), displs(nproc);
  for (int i = 0; i < nproc; i++) {
    displs[i] = distr.range(i).first;
    chunks[i] = distr.range(i).second - distr.range(i).first;
  }
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, first_elem, chunks.data(), displs.data(), MPI_DOUBLE, commun);
#endif
}

} // namespace molpro::linalg::array::util

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_GATHER_ALL_H
