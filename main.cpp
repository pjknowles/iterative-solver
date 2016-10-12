#include "ISDiis.h"
#include "ISDavidson.h"
#include "ISRSPT.h"

extern "C" { void IterativeSolverFTest();}
int main(int argc, char *argv[])
{
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
