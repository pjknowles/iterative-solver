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

  size_t n=10000000;
  size_t repeat=100;
  LinearAlgebra::CachedParameterVector y(n), z(n);
//  IterativeSolver::PagedParameterVector y(n), z(n);
//  IterativeSolver::SimpleParameterVector y(n), z(n);
//  y.setCacheSize(n*1);z.setCacheSize(n*1);
  y.setCacheSize(n-1);z.setCacheSize(n-1);
//  y.setCacheSize(10000);z.setCacheSize(10000);
  y.zero();
  z.zero();
//  IterativeSolver::ParameterScalar one=1;
//  y.put(&one,1,n/2);
//  z.put(&one,1,n/2);
  std::cout <<y.dot(&z)<<std::endl;
  std::clock_t start=std::clock();
//  double result;
  for (size_t r=0; r<repeat; r++) {
//    result += y.dot(&z);
//    y.zero();
    y=z;
//    y.axpy(2.0,&z);
    }
  xout << "time="<<(std::clock()-start)/(double) CLOCKS_PER_SEC<<std::endl;
  if (false) {
//  IterativeSolver::DIIS::randomTest(100,100,0.1,0.0);
//  IterativeSolver::DIIS::randomTest(100,100,0.2,0.0);
//  IterativeSolver::DIIS::randomTest(100,100,0.1,1.0);
//  IterativeSolver::DIIS::randomTest(100,100,0.1,2.0);
  LinearAlgebra::DIIS::randomTest(100,100,0.1,3.0);
//  IterativeSolver::DIIS::test(1,6,1e-10,IterativeSolver::DIIS::KAINmode,0.0002);
  LinearAlgebra::DIIS::test(1,6,1e-10,LinearAlgebra::DIIS::DIISmode,0.2);
  LinearAlgebra::DIIS::test(1,6,1e-10,LinearAlgebra::DIIS::DIISmode,0.0002);
//  IterativeSolver::DIIS::test(1,6,1e6,IterativeSolver::DIIS::disabled,0.0002);
  LinearAlgebra::Davidson::test(9,1,1,2,true);
  LinearAlgebra::Davidson::test(9,1,1,2,false);
  LinearAlgebra::Davidson::test(9,1,1,1);
  LinearAlgebra::Davidson::test(9,1,1,2);
  LinearAlgebra::Davidson::test(9,2,1,2);
  LinearAlgebra::Davidson::test(99,3,1,2);
  LinearAlgebra::RSPT::test(100,2e0);
  IterativeSolverFTest();
    }
#ifdef USE_MPI
  MPI_Finalize();
#endif
  return 0;
}
