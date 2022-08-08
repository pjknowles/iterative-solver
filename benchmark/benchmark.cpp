#include "ArrayBenchmark.h"
#include <iostream>
#ifdef LINEARALGEBRA_ARRAY_GA
#define HAVE_GA_H 1
#endif
#include <molpro/mpi.h>
int main(int argc, char* argv[]) {
  molpro::mpi::init();
#ifdef LINEARALGEBRA_ARRAY_GA
  GA_Initialize();
#endif
  auto rank = molpro::mpi::rank_global();
  auto mpi_size = molpro::mpi::size_global();
  size_t nfast = 100, nslow = 10;
  if (rank == 0)
    std::cout << mpi_size << " MPI ranks" << std::endl;
  for (const auto& length : std::vector<size_t>{500, 1000, 10000, 100000, 1000000, 10000000, 100000000}) {

    if (rank == 0)
      std::cout << "Vector length = " << length << ", numbers of vectors = " << nslow << " / " << nfast;
//    {
//      auto bm = molpro::linalg::ArrayBenchmarkIterable<std::vector<double>>("std::vector<double>", length, nfast, nslow,
//                                                                            false, 0.1);
//      bm.axpy(); // flush memory
//    }
//    if (mpi_size <= 1) {
//      auto bm = molpro::linalg::ArrayBenchmarkIterable<std::vector<double>>("std::vector<double>", length, nfast, nslow,
//                                                                            false, 0.1);
//      bm.all();
//      std::cout << bm;
//    }
#ifdef LINEARALGEBRA_ARRAY_MPI3
    {
      auto bm = molpro::linalg::ArrayBenchmarkDistributed<molpro::linalg::array::DistrArrayMPI3>(
          "DistrArrayMPI3", length, nfast, nslow, false, 0.1);
      if (rank == 0)
        std::cout << ", repeat count = " << bm.m_repeat << std::endl;
      bm.all();
      std::cout << bm;
    }
#else
#endif
#ifdef LINEARALGEBRA_ARRAY_GA
    {
      auto bm = molpro::linalg::ArrayBenchmarkDistributed<molpro::linalg::array::DistrArrayGA>(
          "DistrArrayGA", length, nfast, nslow, false, 0.1);
      bm.all();
      std::cout << bm;
    }
#endif
#ifdef LINEARALGEBRA_ARRAY_HDF5
    {
      auto bm = molpro::linalg::ArrayBenchmarkDDisk<molpro::linalg::array::DistrArrayHDF5>(
          "DistrArrayHDF5 / memory", length, nfast, nslow, false, 0.01);
      bm.all();
      std::cout << bm;
    }
#endif
    {
      auto bm = molpro::linalg::ArrayBenchmarkDDisk<molpro::linalg::array::DistrArrayFile>(
          "DistrArrayFile / memory", length, nfast, nslow, false, 0.1);
      bm.all();
      std::cout << bm;
    }
  }
  molpro::mpi::finalize();
}