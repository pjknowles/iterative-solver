#include "IterativeSolver/ISDiis.h"
#include "IterativeSolver/ISDavidson.h"
#include "IterativeSolver/ISRSPT.h"
#include "IterativeSolver/PagedParameterVector.h"
#include "IterativeSolver/SimpleParameterVector.h"
#include "IterativeSolver/CachedParameterVector.h"
#include <ctime>

extern "C" { void IterativeSolverFTest();}
int main(int argc, char *argv[])
{
  IterativeSolver::CachedParameterVector x(5);
  x.setCacheSize(2);
  for (size_t k=0; k<x.size(); k++) x[k]=k;
  std::cout << "x="<<x.str()<<std::endl;
  std::cout << "x.x="<<x.dot(&x)<<std::endl;

  size_t n=1000000;
  size_t repeat=1000;
  IterativeSolver::CachedParameterVector y(n), z(n);
//  IterativeSolver::PagedParameterVector y(n), z(n);
//  IterativeSolver::SimpleParameterVector y(n), z(n);
  y.zero();
  z.zero();
  y.setCacheSize(n);z.setCacheSize(n);
  y.setCacheSize(1024);z.setCacheSize(1024);
  IterativeSolver::ParameterScalar one=1;
  y.put(&one,1,n/2);
  z.put(&one,1,n/2);
  std::cout <<y.dot(&z)<<std::endl;
  std::clock_t start=std::clock();
  for (size_t r=0; r<repeat; r++)
    y=z;
  xout << "time="<<(std::clock()-start)/(double) CLOCKS_PER_SEC<<std::endl;
  return 0;
//  IterativeSolver::DIIS::randomTest(100,100,0.1,0.0);
//  IterativeSolver::DIIS::randomTest(100,100,0.2,0.0);
//  IterativeSolver::DIIS::randomTest(100,100,0.1,1.0);
//  IterativeSolver::DIIS::randomTest(100,100,0.1,2.0);
  IterativeSolver::DIIS::randomTest(100,100,0.1,3.0);
//  IterativeSolver::DIIS::test(1,6,1e-10,IterativeSolver::DIIS::KAINmode,0.0002);
  IterativeSolver::DIIS::test(1,6,1e-10,IterativeSolver::DIIS::DIISmode,0.2);
  IterativeSolver::DIIS::test(1,6,1e-10,IterativeSolver::DIIS::DIISmode,0.0002);
//  IterativeSolver::DIIS::test(1,6,1e6,IterativeSolver::DIIS::disabled,0.0002);
  IterativeSolver::Davidson::test(9,1,1,2,true);
  IterativeSolver::Davidson::test(9,1,1,2,false);
  IterativeSolver::Davidson::test(9,1,1,1);
  IterativeSolver::Davidson::test(9,1,1,2);
  IterativeSolver::Davidson::test(9,2,1,2);
  IterativeSolver::Davidson::test(99,3,1,2);
  IterativeSolver::RSPT::test(100,2e0);
  IterativeSolverFTest();
  return 0;
}
