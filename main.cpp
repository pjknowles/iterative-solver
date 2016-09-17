#include "Diis.h"
#include "Davidson.h"

int main(int argc, char *argv[])
{
  IterativeSolver::DIIS::test(1);
  IterativeSolver::DIIS::test(1,6,1e6,IterativeSolver::DIIS::DIISmode,0.2);
//  IterativeSolver::DIIS::test(1,6,1e6,IterativeSolver::DIIS::disabled,0.0002);
  IterativeSolver::DIIS::test(1,6,1e6,IterativeSolver::DIIS::DIISmode,0.4);
  IterativeSolver::DIIS::test(1,6,1e6,IterativeSolver::DIIS::DIISmode,0.8);
  IterativeSolver::DIIS::test(1,6,1e6,IterativeSolver::DIIS::DIISmode,0.8);
//  IterativeSolver::DIIS::test(1,6,1e6,IterativeSolver::DIIS::KAINmode);
  return 0;
}
