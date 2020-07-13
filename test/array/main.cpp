#ifdef TEST_WITH_GA
#include <ga-mpi.h>
#include <ga.h>
#include <macdecls.h>
#endif
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <mpi.h>

#include "parallel_util.h"

MPI_Comm molpro::linalg::test::mpi_comm;

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  MPI_Init(&argc, &argv);
  molpro::linalg::test::mpi_comm = MPI_COMM_WORLD;
#ifdef TEST_WITH_GA
  GA_Initialize();
  int mem = 10000000;
  MA_init(C_DBL, mem, mem);
  if (!GA_Create_mutexes(1))
    GA_Error((char *)"Failed to create mutexes", 1);
  molpro::linalg::test::mpi_comm = GA_MPI_Comm();
#endif
  int result = RUN_ALL_TESTS();
#ifdef TEST_WITH_GA
  GA_Destroy_mutexes();
  GA_Terminate();
#endif
  MPI_Finalize();
  return result;
}
