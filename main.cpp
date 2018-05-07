#include "IterativeSolver/ISDiis.h"
#include "IterativeSolver/ISDavidson.h"
#include "IterativeSolver/ISRSPT.h"
#include "IterativeSolver/PagedParameterVector.h"
#include "IterativeSolver/SimpleParameterVector.h"
#include "IterativeSolver/CachedParameterVector.h"
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
    if (rank==0) std::cout << size<<" MPI processs"<<(size>1?"es ":"")<<std::endl;
//    std::cout << "MPI process number "<<rank<<std::endl;
  }
#endif
//  IterativeSolver::CachedParameterVector x(5);
//  x.setCacheSize(2);
//  for (size_t k=0; k<x.size(); k++) x[k]=k;
//  std::cout << "x="<<x.str()<<std::endl;
//  std::cout << "x.x="<<x.dot(&x)<<std::endl;

  if (true) {
//  IterativeSolver::DIIS::randomTest(100,100,0.1,0.0);
//  IterativeSolver::DIIS::randomTest(100,100,0.2,0.0);
//  IterativeSolver::DIIS::randomTest(100,100,0.1,1.0);
//  IterativeSolver::DIIS::randomTest(100,100,0.1,2.0);
//  LinearAlgebra::DIIS<double>::randomTest(100,100,0.1,3.0);
//  IterativeSolver::DIIS::test(1,6,1e-10,IterativeSolver::DIIS::KAINmode,0.0002);
//  LinearAlgebra::DIIS<double>::test(1,6,1e-10,LinearAlgebra::DIIS<double>::DIISmode,0.2);
//  LinearAlgebra::DIIS<double>::test(1,6,1e-10,LinearAlgebra::DIIS<double>::DIISmode,0.0002);
//  IterativeSolver::DIIS::test(1,6,1e6,IterativeSolver::DIIS::disabled,0.0002);
  LinearAlgebra::Davidson<double>::test(2,1,1,2,true);
  LinearAlgebra::Davidson<double>::test(9,1,1,2,true);
  LinearAlgebra::Davidson<double>::test(9,1,1,2,false);
  LinearAlgebra::Davidson<double>::test(9,1,1,1);
  LinearAlgebra::Davidson<double>::test(9,1,1,2);
  LinearAlgebra::Davidson<double>::test(9,2,1,2);
  LinearAlgebra::Davidson<double>::test(99,3,1,2);
//  LinearAlgebra::RSPT<double>::test(100,2e0);
  IterativeSolverFTest();
    }
#ifdef USE_MPI
  MPI_Finalize();
#endif
  return 0;
}
