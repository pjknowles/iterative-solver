#include "ArrayBenchmark.h"
#include <iostream>
#include <molpro/mpi.h>
#ifdef HAVE_MPI_H
#include <mpi.h>
#endif
int main(int argc, char* argv[]) {
  molpro::mpi::init();
  auto rank = molpro::mpi::rank_global();
  auto mpi_size = molpro::mpi::size_global();
//#ifdef LINEARALGEBRA_ARRAY_GA
//  GA_Initialize();
//#endif
  if (rank == 0)
    std::cout << mpi_size << " MPI ranks" << std::endl;
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
  molpro::mpi::finalize();
}