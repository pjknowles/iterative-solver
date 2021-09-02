#ifndef LINEARALGEBRA_TEST_ARRAY_PARALLEL_UTIL_H
#define LINEARALGEBRA_TEST_ARRAY_PARALLEL_UTIL_H
#include <mpi.h>
#include <numeric>

#include <molpro/linalg/array/DistrArrayFile.h>
#include <molpro/linalg/array/DistrArraySpan.h>
#include <molpro/linalg/array/Span.h>
#include <molpro/linalg/array/util/Distribution.h>

namespace molpro {
namespace linalg {
namespace test {

extern MPI_Comm mpi_comm;

std::tuple<std::vector<molpro::linalg::array::DistrArraySpan>, std::vector<molpro::linalg::array::DistrArraySpan>,
           std::vector<molpro::linalg::array::DistrArrayFile>>
get_contiguous(size_t n, size_t dim);

} // namespace test
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_TEST_ARRAY_PARALLEL_UTIL_H
