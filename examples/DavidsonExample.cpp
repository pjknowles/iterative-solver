#include "ISDavidson.h"
#include "SimpleParameterVector.h"

using namespace IterativeSolver;

// Find lowestG eigensolutions of M(i,j) = alpha*(i+1)*delta(i,j) + i + j
static double n; // dimension of problem
static double alpha; // separation of diagonal elements

static void _matrix_residual(const ParameterVectorSet & psx, ParameterVectorSet & outputs, std::vector<ParameterScalar> shift=std::vector<ParameterScalar>(), bool append=false) {
  for (size_t k=0; k<psx.size(); k++) {
      for (size_t i=0; i<n; i++) {
          if (append)
            (*outputs[k])[i] += (alpha*(i+1))*(*psx[k])[i];
          else
            (*outputs[k])[i] = (alpha*(i+1))*(*psx[k])[i];
          for (size_t j=0; j<n; j++)
            (*outputs[k])[i] += (i+j)*(*psx[k])[j];
        }
    }
}

static void _matrix_preconditioner(const ParameterVectorSet & psg, ParameterVectorSet & psc, std::vector<ParameterScalar> shift=std::vector<ParameterScalar>(), bool append=false) {
  for (size_t k=0; k<psc.size(); k++) {
      if (append) {
          for (size_t i=0; i<n; i++)
            (*psc[k])[i] -= (*psg[k])[i]/(shift[k]+2*i+alpha*(i+1));
        } else {
          for (size_t i=0; i<n; i++)
            (*psc[k])[i] =- (*psg[k])[i]/(shift[k]+2*i+alpha*(i+1));
        }
    }
}

int main(int argc, char *argv[])
{
  alpha=100;
  n=10;
  Davidson solver(&_matrix_residual,&_matrix_preconditioner);
  solver.m_verbosity=1;
  solver.m_roots=2;
  ParameterVectorSet g;
  ParameterVectorSet x;
  for (size_t root=0; root<solver.m_roots; root++) {
      SimpleParameterVector* xx = new SimpleParameterVector(n); x.push_back(xx);
      SimpleParameterVector* gg = new SimpleParameterVector(n); g.push_back(gg);
      xx->zero(); (*xx)[root]=1; // initial guess
    }
  if (not solver.solve(g,x)) std::cout << "Failure"<<std::endl;
  std::cout << "Error=";
  for (size_t root=0; root<solver.m_roots; root++)
    std::cout <<solver.errors()[root]<<" ";
  std::cout <<"after "<<solver.iterations()<<" iterations"<<std::endl;
  for (size_t root=0; root<solver.m_roots; root++) {
      std::cout << "Eigenvector:"; for (size_t k=0; k<n; k++) std::cout<<" "<<(*x[root])[k]; std::cout<<std::endl;
      delete x[root];
      delete g[root];
    }
}
