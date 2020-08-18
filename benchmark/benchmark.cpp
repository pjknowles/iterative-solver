#include "ArrayBenchmark.h"
#ifdef HAVE_MPI_H
#include <mpi.h>
#endif
int main(int argc, char* argv[]) {
  int rank = 0;
#ifdef HAVE_MPI_H
  MPI_Init(&argc, &argv);
  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0)
    std::cout << mpi_size << " MPI ranks" << std::endl;
#endif
#ifdef LINEARALGEBRA_ARRAY_GA
  GA_Initialize();
#endif
  for (const auto& length : std::vector<size_t>{100000, 1000000, 10000000, 100000000}) {

    std::cout << "Vector length " << length << std::endl;
    molpro::linalg::ArrayBenchmarkII<std::vector<double>>("std::vector<double>", length, true, 1e-6)
        .axpy(); // flush memory
    if (mpi_size <= 1) {
      auto bm = molpro::linalg::ArrayBenchmarkII<std::vector<double>>("std::vector<double>", length);
      bm.all();
      if (rank == 0)
        std::cout << bm;
    }
#ifdef LINEARALGEBRA_ARRAY_MPI3
    {
      auto bm = molpro::linalg::ArrayBenchmarkDD<molpro::linalg::array::DistrArrayMPI3>(
          "molpro::linalg::array::DistrArrayMPI3", length);
      bm.all();
      std::cout << bm;
    }
//    {
//      auto bm = molpro::linalg::ArrayBenchmarkID<std::vector<double>, molpro::linalg::array::DistrArrayMPI3>(
//          "std::vector<double>,molpro::linalg::array::DistrArrayMPI3", length);
//      bm.all();
//      std::cout << bm;
//    }
#endif
#ifdef LINEARALGEBRA_ARRAY_GA
    {
      auto bm = molpro::linalg::ArrayBenchmarkDD<molpro::linalg::array::DistrArrayGA>(
          "molpro::linalg::array::DistrArrayGA", length);
      bm.all();
      std::cout << bm;
    }
#endif
#ifdef LINEARALGEBRA_ARRAY_HDF5
    {
      //    auto bm = molpro::linalg::ArrayBenchmarkDD<molpro::linalg::array::DistrArrayHDF5>(
      //        "molpro::linalg::array::DistrArrayHDF5",length);
      //    bm.all();
      //    std::cout << bm ;
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