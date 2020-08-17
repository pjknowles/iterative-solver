#include "ArrayBenchmark.h"
#ifdef HAVE_MPI_H
#include <mpi.h>
#endif
int main(int argc, char* argv[]) {
#ifdef HAVE_MPI_H
  MPI_Init(&argc, &argv);
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  std::cout << size << std::endl;
#endif
#ifdef LINEARALGEBRA_ARRAY_GA
  GA_Initialize();
#endif
  {

    auto bm = molpro::linalg::ArrayBenchmarkII<std::vector<double>>("std::vector<double>");
    bm.all();
    std::cout << bm << std::endl;
  }
#ifdef LINEARALGEBRA_ARRAY_MPI3
  {
    auto bm = molpro::linalg::ArrayBenchmarkDD<molpro::linalg::array::DistrArrayMPI3>(
        "molpro::linalg::array::DistrArrayMPI3");
    bm.all();
    std::cout << bm << std::endl;
  }
#endif
#ifdef LINEARALGEBRA_ARRAY_GA
  {
    auto bm =
        molpro::linalg::ArrayBenchmarkDD<molpro::linalg::array::DistrArrayGA>("molpro::linalg::array::DistrArrayGA");
    bm.all();
    std::cout << bm << std::endl;
  }
#endif
#ifdef LINEARALGEBRA_ARRAY_HDF5
  {
    //    auto bm = molpro::linalg::ArrayBenchmarkDD<molpro::linalg::array::DistrArrayHDF5>(
    //        "molpro::linalg::array::DistrArrayHDF5");
    //    bm.all();
    //    std::cout << bm << std::endl;
  }
#endif
#ifdef HAVE_MPI_H
  MPI_Finalize();
#endif
}