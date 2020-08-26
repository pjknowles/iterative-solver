#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_GATHER_ALL_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_GATHER_ALL_H
#include <mpi.h>
#include "Distribution.h"

namespace molpro {
namespace linalg {
namespace array {
namespace util {

/*!
 * @brief Replicate data of a full container on all the processes based on distributed pieces
 *
 * @tparam X type of the container, must have range() member function, specifying local begin and end indices
 * @param commun MPI communicator
 * @param x container, which data to be replicated
 */
template<class X>
void gather_all(X& x, MPI_Comm commun) {
  int nproc;
  MPI_Comm_size(commun, &nproc);
  int chunks[nproc], displs[nproc];
  for (int i = 0; i < nproc; i++) {
    displs[i] = x.distribution().range(i).first;
    chunks[i] = x.distribution().range(i).second - x.distribution().range(i).first;;
  }
  if (std::is_same<typename X::value_type, double>::value) {
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &(*x.local_buffer())[0], chunks, displs, MPI_DOUBLE, commun);
  } else {
    // TODO: through an error? case for std::complex?
  }
}

} // namespace util
} // namespace array
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_GATHER_ALL_H
