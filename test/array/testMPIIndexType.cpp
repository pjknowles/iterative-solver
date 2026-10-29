#include <gtest/gtest.h>

#include <molpro/linalg/array/DistrArray.h>

#include <mpi.h>

using molpro::linalg::array::DistrArray;
using molpro::linalg::array::util::mpi_index_type;

TEST(mpi_index_type, size_matches_index_type) {
  int size = 0;
  MPI_Type_size(mpi_index_type(), &size);
  EXPECT_EQ(static_cast<size_t>(size), sizeof(DistrArray::index_type));
}
