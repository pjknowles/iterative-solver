#include "ArrayBenchmark.h"
#ifdef HAVE_MPI_H
#include <mpi.h>
#endif
int main(int argc, char* argv[]) {
  int rank = 0;
#ifdef HAVE_MPI_H
  MPI_Init(&argc, &argv);
#ifdef LINEARALGEBRA_ARRAY_GA
  GA_Initialize();
#endif
  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0)
    std::cout << mpi_size << " MPI ranks" << std::endl;
#else
  int mpi_size = 1;
#endif
  for (const auto& length : std::vector<size_t>{500, 1000, 10000, 100000, 1000000, 10000000, 100000000}) {

    if (rank == 0)
      std::cout << "Vector length = " << length;
    {
      auto bm = molpro::linalg::ArrayBenchmarkIterable<std::vector<double>>("std::vector<double>", length, false, 0.1);
      bm.axpy(); // flush memory
      if (rank == 0)
        std::cout << ", repeat count = " << bm.m_repeat << std::endl;
    }
    if (mpi_size <= 1) {
      auto bm = molpro::linalg::ArrayBenchmarkIterable<std::vector<double>>("std::vector<double>", length, false, 0.1);
      bm.all();
      std::cout << bm;
    }
#ifdef LINEARALGEBRA_ARRAY_MPI3
    {
      auto bm = molpro::linalg::ArrayBenchmarkDistributed<molpro::linalg::array::DistrArrayMPI3>(
          "molpro::linalg::array::DistrArrayMPI3", length, false, 0.1);
      bm.all();
      std::cout << bm;
    }
#endif
#ifdef LINEARALGEBRA_ARRAY_GA
    {
      auto bm = molpro::linalg::ArrayBenchmarkDistributed<molpro::linalg::array::DistrArrayGA>(
          "molpro::linalg::array::DistrArrayGA", length, false, 0.1);
      bm.all();
      std::cout << bm;
    }
#endif
#ifdef LINEARALGEBRA_ARRAY_HDF5
    {
      auto bm = molpro::linalg::ArrayBenchmarkDDisk<molpro::linalg::array::DistrArrayHDF5>(
          "molpro::linalg::array::DistrArrayHDF5", length, false, 0.01);
      bm.all();
      std::cout << bm;
    }
#endif
  }
#ifdef LINEARALGEBRA_ARRAY_GA
  GA_Terminate();
#endif
#ifdef HAVE_MPI_H
  MPI_Finalize();
#endif
}