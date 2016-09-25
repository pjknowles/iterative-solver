#include "Diis.h"
#include "Davidson.h"

extern "C" { void IterativeSolverFTest();}
int main(int argc, char *argv[])
{
  IterativeSolver::DIIS::test(1);
  IterativeSolver::DIIS::test(1,6,1e6,IterativeSolver::DIIS::DIISmode,0.2);
//  IterativeSolver::DIIS::test(1,6,1e6,IterativeSolver::DIIS::disabled,0.0002);
  IterativeSolver::DIIS::test(1,6,1e6,IterativeSolver::DIIS::DIISmode,0.4);
  IterativeSolver::DIIS::test(1,6,1e6,IterativeSolver::DIIS::DIISmode,0.8);
  IterativeSolver::DIIS::test(1,6,1e6,IterativeSolver::DIIS::DIISmode,0.8);
//  IterativeSolver::DIIS::test(1,6,1e6,IterativeSolver::DIIS::KAINmode);
  IterativeSolver::Davidson::test(9,1,1,0);
  IterativeSolver::Davidson::test(9,1,1,1);
  IterativeSolver::Davidson::test(9,1,1,2);
  IterativeSolver::Davidson::test(9,2,1,2);
  IterativeSolver::Davidson::test(99,3,1,2);
  IterativeSolverFTest();
  return 0;
}
