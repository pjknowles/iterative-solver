#include "IterativeSolver/ISDiis.h"
#include "IterativeSolver/SimpleParameterVector.h"

using namespace IterativeSolver;

static double alpha;
static double anharmonicity;
static double n;

static void _anharmonic_residual(const ParameterVectorSet & psx, ParameterVectorSet & outputs, std::vector<ParameterScalar> shift=std::vector<ParameterScalar>(), bool append=false) {
  std::vector<ParameterScalar> psxk(n);
  std::vector<ParameterScalar> output(n);
  psx.front()->get(&(psxk[0]),n,0);
  if (append)
    outputs.front()->get(&(output[0]),n,0);
  else
    outputs.front()->zero();
  for (size_t i=0; i<n; i++) {
      output[i] += (alpha*(i+1)+anharmonicity*psxk[i])*psxk[i];
      for (size_t j=0; j<n; j++)
        output[i] += (i+j)*psxk[j];
    }
  outputs.front()->put(&output[0],n,0);
}
static void _anharmonic_preconditioner(const ParameterVectorSet & psg, ParameterVectorSet & psc, std::vector<ParameterScalar> shift=std::vector<ParameterScalar>(), bool append=false) {
  std::vector<ParameterScalar> psck(n);
  std::vector<ParameterScalar> psgk(n);
  psg.front()->get(&psgk[0],n,0);
  if (append) {
      psc.front()->get(&psck[0],n,0);
      for (size_t i=0; i<n; i++)
        psck[i] -= psgk[i]/(alpha*(i+1));
    } else {
      for (size_t i=0; i<n; i++)
        psck[i] =- psgk[i]/(alpha*(i+1));
    }
  psc.front()->put(&psck[0],n,0);
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
  xx.zero(); double one=1; xx.put(&one,1,0);  // initial guess
  if (not solver.solve(g,x)) std::cout << "Failure"<<std::endl;
  std::cout << "Distance of solution from origin: "<<std::sqrt(xx.dot(&xx))<<std::endl;
  std::cout << "Error="<<solver.errors().front()<<" after "<<solver.iterations()<<" iterations"<<std::endl;
}
