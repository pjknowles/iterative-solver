#include "IterativeSolver/ISDiis.h"
#include "IterativeSolver/ISDavidson.h"
#include "IterativeSolver/ISRSPT.h"
#include "PagedVector.h"
#include <ctime>
#ifdef USE_MPI
#include <mpi.h>
#endif

extern "C" { void IterativeSolverFTest();}
int main(int argc, char *argv[])
{
#ifdef USE_MPI
  MPI_Init(&argc,&argv);
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    if (rank==0) std::cout << size<<" MPI process"<<(size>1?"es ":"")<<std::endl;
//    std::cout << "MPI process number "<<rank<<std::endl;
  }
#endif

  if (true) {
   using namespace LinearAlgebra;
//  IterativeSolver::DIIS::randomTest(100,100,0.1,0.0);
//  IterativeSolver::DIIS::randomTest(100,100,0.2,0.0);
//  IterativeSolver::DIIS::randomTest(100,100,0.1,1.0);
//  IterativeSolver::DIIS::randomTest(100,100,0.1,2.0);
//  LinearAlgebra::DIIS<double>::randomTest(100,100,0.1,3.0);
//  DIISTest<PagedVector<double> >(1,6,1e-10,LinearAlgebra::DIIS<double>::DIISmode,0.0002);
//  DIISTest<PagedVector<double> >(1,6,1e-10,LinearAlgebra::DIIS<double>::DIISmode,0.2);
//  DIISTest<PagedVector<double> >(1,6,1e-3,LinearAlgebra::DIIS<double>::disabled,0.0002);
  DavidsonTest<PagedVector<double> >(2,1,1,2,true);
  DavidsonTest<PagedVector<double> >(9,1,1,2,true);
  DavidsonTest<PagedVector<double> >(9,1,1,2,false);
  DavidsonTest<PagedVector<double> >(9,1,1,1);
  DavidsonTest<PagedVector<double> >(9,1,1,2);
  DavidsonTest<PagedVector<double> >(9,2,1,2);
  DavidsonTest<PagedVector<double> >(100,1,1,2);
  DavidsonTest<PagedVector<double> >(100,3,1,2);
  DavidsonTest<PagedVector<double> >(100,3,1,2,false);
  DavidsonTest<PagedVector<double> >(600,3,1,2,true);
//  RSPTTest<PagedVector<double> ,double>(100,2e0);
  IterativeSolverFTest();
    }
#ifdef USE_MPI
  MPI_Finalize();
#endif
  return 0;
}
