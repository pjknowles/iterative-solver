#include <iostream>
#include <vector>
#include <cmath>
#include "Diis.h"

using namespace std;
using namespace IterativeSolver;

typedef std::vector<double> doubles;

bool rosenbrock=true;


void residual(const ParameterVectorSet & psx, ParameterVectorSet & outputs, std::vector<ParameterScalar> shift=std::vector<ParameterScalar>()) {
    ParameterVector x=psx.front();
    std::cout << "in residual psx="<<psx<<std::endl;
    std::cout << "in residual x="<<x<<std::endl;
    ParameterVector result(2);
	  std::cout << "in residual x="<<x[0]<<","<<x[1]<<std::endl;
	  if (rosenbrock) {
      result[0]=(2*x[0]-2+400*x[0]*(x[0]*x[0]-x[1])); result[1]=(200*(x[1]-x[0]*x[0])); // Rosenbrock
	    } else {
      result[0]=(x[0]-1);
      result[1]=(x[1]-1);
	    }
	  std::cout << "in residual result="<<result[0]<<","<<result[1]<<" size="<<result.size()<<std::endl;
      outputs.front()=result;
}

void updater(const ParameterVectorSet & psg, ParameterVectorSet & psc, std::vector<ParameterScalar> shift) {
       ParameterVector c=psc.front();
       ParameterVector g=psg.front();
    ParameterVector diag(2);
    diag[0]=700; diag[1]=200;
       std::cout << "updater receives "<<psc<<psg<<std::endl;
    if (rosenbrock) {
        if (false) {
        std::cout << "preconditioner "<<1/(2+1200*c[0]*c[0]-400*c[1])<<", "<<1.0/200<<std::endl;
        c[0]-=g[0]/(2+1200*c[0]*c[0]-400*c[1]);
		c[1]-=g[1]/200;
        } else
        {
            std::cout << "preconditioner "<<1/diag[0]<<", "<<1.0/diag[1]<<std::endl;
            c[0] -=g[0]/diag[0];
            c[1] -=g[1]/diag[1];
        }
    } else {
        c[0]-=g[0];
		c[1]-=g[1];
	}
    psc.front()=c;
       std::cout << "updater returns "<<psc<<std::endl;
}

void simple() {
  ParameterVectorSet x; x.push_back(ParameterVector(2));
  x.front()[0]=x.front()[1]=0.9; // initial guess
  std::cout << "after constructing x "<<x<<std::endl;
  ParameterVectorSet g; g.push_back(ParameterVector(2));
  std::cout << "after constructing g "<<g<<std::endl;
  Diis d(&updater,&residual);
  bool converged=false;
  for (int iteration=1; iteration < 15 && not converged; iteration++) {
  std::cout << "before residual, x "<<x<<std::endl;
      residual(x,g);
      converged = d.iterate(g,x);
      std::cout << "iteration "<<iteration<<", Residual norm = "<<std::sqrt(d.fLastResidual())
                  << ", Distance from solution = "<<std::sqrt((x.front()[0]-1)*(x.front()[0]-1)+(x.front()[1]-1)*(x.front()[1]-1))
              <<", converged? "<<converged
      <<std::endl;
    }
}

int main(int argc, char *argv[])
{
  cout << "Test DIIS" << endl;
  std::cout << "try simple version"<<std::endl;
  simple();
  return 0;
}
