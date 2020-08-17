#include "ArrayBenchmark.h"
#ifdef HAVE_MPI_H
#include <mpi.h>
#endif
int main(int argc, char* argv[]) {
#ifdef HAVE_MPI_H
  MPI_Init(&argc, &argv);
#endif
  molpro::linalg::ArrayBenchmark<std::vector<double>> bm("std::vector<double>");
  bm.all();
  std::cout << bm << std::endl;
#ifdef HAVE_MPI_H
  MPI_Finalize();
#endif
}