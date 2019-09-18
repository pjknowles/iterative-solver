#define HAVE_MPI_H
#include <iostream>
#include "PagedVector.h"
#include <chrono>
#include <mpi.h>
constexpr size_t window_size = 10000;
constexpr size_t vector_size = 50000000;
int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  int mpiSize;
  MPI_Comm_size(MPI_COMM_COMPUTE, &mpiSize);
  int mpiRank;
  MPI_Comm_rank(MPI_COMM_COMPUTE, &mpiRank);
  if (mpiRank==0) std::cout << mpiSize<<" MPI processes"<<std::endl;
  using pv = LinearAlgebra::PagedVector<double, window_size>;
  for (int option = 0; option < 4; option++) {
    pv v0(vector_size,0);
    MPI_Barrier(MPI_COMM_COMPUTE);
    auto startTiming = std::chrono::steady_clock::now();
    auto v1 = pv(v0, option);
    auto endTiming = std::chrono::steady_clock::now();
    MPI_Barrier(MPI_COMM_COMPUTE);
    if (mpiRank==0)
    std::cout << " copy PagedVector from memory, construction option="<<option<<":  seconds="
              << std::chrono::duration_cast<std::chrono::nanoseconds>(endTiming - startTiming).count() * 1e-9
              << ", bandwidth="
        << vector_size*sizeof(double)/double(std::chrono::duration_cast<std::chrono::nanoseconds>(endTiming - startTiming).count())<<" Gb/s"
              << std::endl;
    if (v0.dot(v0) != v1.dot(v1)) throw std::runtime_error("copy error");
    if (v0.dot(v1) != v1.dot(v1)) throw std::runtime_error("copy error");
  }
  MPI_Finalize();
}
