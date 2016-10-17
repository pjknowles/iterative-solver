#include "ISDiis.h"
#include "SimpleParameterVector.h"

using namespace IterativeSolver;

static double alpha;
static double anharmonicity;
static double n;

static void _anharmonic_residual(const ParameterVectorSet & psx, ParameterVectorSet & outputs, std::vector<ParameterScalar> shift=std::vector<ParameterScalar>(), bool append=false) {
  if (not append) outputs.front()->zero();
  for (size_t i=0; i<n; i++) {
      (*outputs.front())[i] += (alpha*(i+1)+anharmonicity*(*psx.front())[i])*(*psx.front())[i];
      for (size_t j=0; j<n; j++)
        (*outputs.front())[i] += (i+j)*(*psx.front())[j];
    }
}
static void _anharmonic_preconditioner(const ParameterVectorSet & psg, ParameterVectorSet & psc, std::vector<ParameterScalar> shift=std::vector<ParameterScalar>(), bool append=false) {
  if (append) {
      for (size_t i=0; i<n; i++)
        (*psc.front())[i] -= (*psg.front())[i]/(alpha*(i+1));
    } else {
      for (size_t i=0; i<n; i++)
        (*psc.front())[i] =- (*psg.front())[i]/(alpha*(i+1));
    }
}

int main(int argc, char *argv[])
{
  alpha=1;
  n=100;
  anharmonicity=.5;
  DIIS solver(&_anharmonic_residual,&_anharmonic_preconditioner);
  solver.m_verbosity=1;
  SimpleParameterVector gg(n); ParameterVectorSet g; g.push_back(&gg);
  SimpleParameterVector xx(n); ParameterVectorSet x; x.push_back(&xx);
  xx.zero(); xx[n]=xx[0]=1; // initial guess
  if (not solver.solve(g,x)) std::cout << "Failure"<<std::endl;
  std::cout << "Distance of solution from origin: "<<std::sqrt(xx.dot(&xx))<<std::endl;
  std::cout << "Error="<<solver.errors().front()<<" after "<<solver.iterations()<<" iterations"<<std::endl;
}
