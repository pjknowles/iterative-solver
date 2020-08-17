#include "ArrayBenchmark.h"
#ifdef HAVE_MPI_H
#include <mpi.h>
#endif
int main(int argc, char* argv[]) {
#ifdef HAVE_MPI_H
  MPI_Init(&argc, &argv);
#endif
  {

    auto bm = molpro::linalg::ArrayBenchmarkII<std::vector<double>>("std::vector<double>");
    bm.all();
    std::cout << bm << std::endl;
  }
#ifdef HAVE_MPI_H
  {
    auto bm = molpro::linalg::ArrayBenchmarkDD<molpro::linalg::array::DistrArrayMPI3>(
        "molpro::linalg::array::DistrArrayMPI3");
    bm.all();
    std::cout << bm << std::endl;
  }
  MPI_Finalize();
#endif
}