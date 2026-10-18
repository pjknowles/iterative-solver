#include "ArrayBenchmark.h"
#include <iostream>
#ifdef LINEARALGEBRA_ARRAY_GA
#define HAVE_GA_H 1
#endif
#include <molpro/linalg/options.h>
#include <molpro/mpi.h>
int main(int argc, char* argv[]) {
  molpro::mpi::init();
#ifdef LINEARALGEBRA_ARRAY_GA
  GA_Initialize();
#endif
  auto rank = molpro::mpi::rank_global();
  auto mpi_size = molpro::mpi::size_global();
  const size_t nSlow = 100, nFast = 10;
  molpro::linalg::set_options(molpro::Options("ITERATIVE-SOLVER", "GEMM_PAGESIZE=8192, GEMM_BUFFERS=2" ));
  std::cout << "GEMM_PAGESIZE=" << molpro::linalg::options()->parameter("GEMM_PAGESIZE", 0) << std::endl;
  std::cout << "GEMM_BUFFERS=" << molpro::linalg::options()->parameter("GEMM_BUFFERS", 0) << std::endl;
  if (rank == 0)
    std::cout << mpi_size << " MPI ranks" << std::endl;
  for (const auto& length : std::vector<size_t>{500, 1000, 10000, 100000, 1000000}) {

    if (rank == 0)
      std::cout << "Vector length = " << length << ", numbers of vectors = " << nFast << " / " << nSlow << std::endl;

#ifdef LINEARALGEBRA_ARRAY_MPI3
    {
      auto bm = molpro::linalg::ArrayBenchmarkDistributed<molpro::linalg::array::DistrArrayMPI3>(
          "DistrArrayMPI3", length, nSlow, nFast, false, 0.1);
      bm.all();
      std::cout << bm;
      bm.profiler().dotgraph(bm.m_title + "." + std::to_string(length) + ".gv");
    }
#else
#endif
#ifdef LINEARALGEBRA_ARRAY_GA
    // TODO: It's presently not possible to do this because DistrArrayGA needs to know the communicator?
    // TODO: consider writing a resource class to allow storage of information needed for construction of any DistrArray
//    {
//      auto bm = molpro::linalg::ArrayBenchmarkDistributed<molpro::linalg::array::DistrArrayGA>(
//          "DistrArrayGA", length, nSlow, nFast, false, 0.1);
//      bm.all();
//      std::cout << bm;
//    }
#endif

    if (true) {
      auto bm = molpro::linalg::ArrayBenchmarkDDisk<molpro::linalg::array::DistrArrayFile>("DistrArrayFile", length,
                                                                                           nSlow, nFast, false, 0.1);
      bm.all();
      std::cout << bm;
      bm.profiler().dotgraph(bm.m_title + "." + std::to_string(length) + ".gv");
    }
#ifdef LINEARALGEBRA_ARRAY_HDF5
    if (true) {
      auto bm = molpro::linalg::ArrayBenchmarkDDisk<molpro::linalg::array::DistrArrayHDF5>("DistrArrayHDF5", length,
                                                                                           nSlow, nFast, false, 0.01);
      bm.all();
      std::cout << bm;
      bm.profiler().dotgraph(bm.m_title + "." + std::to_string(length) + ".gv");
    }
#endif
  }

  // Operations in both directions between arrays in memory and on disk: std::vector and Span with the handlers for
  // contiguous arrays, compared with DistrArraySpan with the handlers for distributed arrays
  for (const auto& length : std::vector<size_t>{500, 1000, 10000, 100000, 1000000}) {
    if (rank == 0)
      std::cout << "Vector length = " << length << ", numbers of vectors = " << nFast << " / " << nSlow << std::endl;
    {
      auto bm = molpro::linalg::ArrayBenchmarkDDisk<molpro::linalg::array::DistrArrayFile>(
          "DistrArraySpan-DistrArrayFile", length, nSlow, nFast, false, 0.1);
      bm.all_both_directions();
      std::cout << bm;
    }
    {
      auto bm = molpro::linalg::ArrayBenchmarkIterableDDisk<molpro::linalg::array::DistrArrayFile>(
          "vector-DistrArrayFile", length, nSlow, nFast, false, 0.1);
      bm.all_both_directions();
      std::cout << bm;
    }
    {
      // views into one buffer, so evenly spaced, unlike separate std::vector
      auto bm = molpro::linalg::ArrayBenchmarkIterableDDisk<molpro::linalg::array::DistrArrayFile,
                                                            molpro::linalg::array::Span<double>>(
          "Span-DistrArrayFile", length, nSlow, nFast, false, 0.1);
      bm.all_both_directions();
      std::cout << bm;
    }
  }

  // Compare accessing a DistrArray through the generic STL iterator interface against its own
  // native, optimized interface (local_buffer()/fill()). Kept to modest lengths since the
  // iterator does one remote-memory-access call per element.
  for (const auto& length : std::vector<size_t>{100, 1000, 10000}) {
    if (rank == 0)
      std::cout << "STL iterator vs. native access, vector length = " << length << std::endl;
#ifdef LINEARALGEBRA_ARRAY_MPI3
    molpro::linalg::benchmarkIteratorVsNative<molpro::linalg::array::DistrArrayMPI3>(
        "DistrArrayMPI3.iterator_vs_native." + std::to_string(length), length, 0.2);
#endif
#ifdef LINEARALGEBRA_ARRAY_GA
    molpro::linalg::benchmarkIteratorVsNative<molpro::linalg::array::DistrArrayGA>(
        "DistrArrayGA.iterator_vs_native." + std::to_string(length), length, 0.2);
#endif
  }
#ifdef LINEARALGEBRA_ARRAY_GA
  GA_Terminate();
#endif
  molpro::mpi::finalize();
}