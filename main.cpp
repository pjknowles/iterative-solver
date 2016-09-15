#include <iostream>
#include <vector>
#include <cmath>
#include "Diis.h"

using namespace std;
using namespace IterativeSolver;

typedef std::vector<double> doubles;

bool rosenbrock=true;


void residual(const ParameterVectorSet & psx, ParameterVectorSet & outputs, std::vector<ParameterScalar> shift=std::vector<ParameterScalar>()) {
	  doubles x(2); ParameterVector px(&x[0],x.size()); px=psx.front();
	  doubles result;
	  std::cout << "in residual x="<<x[0]<<","<<x[1]<<std::endl;
	  if (rosenbrock) {
	  result.push_back(2*x[0]-2+400*x[0]*(x[0]*x[0]-x[1])); result.push_back(200*(x[1]-x[0]*x[0])); // Rosenbrock
	    } else {
	  result.push_back(x[0]-1);
	  result.push_back(x[1]-1);
	    }
	  std::cout << "in residual result="<<result[0]<<","<<result[1]<<" size="<<result.size()<<std::endl;
	  ParameterVector pr(&result[0],result.size());
	  std::cout << "pr: "<<pr<<" size="<<pr.size()<<std::endl;
	  outputs.front()=pr;
}

void updater(const ParameterVectorSet & psg, ParameterVectorSet & psc, std::vector<ParameterScalar> shift) {
	  doubles c(2); ParameterVector pc(&c[0],c.size()); pc=psc.front();
	  doubles g(2); ParameterVector pg(&g[0],g.size()); pg=psg.front();
	if (rosenbrock) {
		c[0]-=g[0]/(2+1200*c[0]*c[0]-400*c[1]);
		c[1]-=g[1]/200;
	} else {
		c[0]-=g[0];
		c[1]-=g[1];
	}
	psc.front()=pc;
}

void simple() {
  std::vector<double> xx(2);
  std::vector<double> gg(2);
  xx[0]=xx[1]=0.9; // initial guess
  std::cout << "before constructing x"<<std::endl;
  ParameterVectorSet x; x.push_back(ParameterVector(&xx[0],xx.size()));
  std::cout << "after constructing x"<<std::endl;
  ParameterVectorSet g; g.push_back(ParameterVector(&gg[0],gg.size()));
  std::cout << "after constructing g"<<std::endl;
  Diis d(&updater,&residual);
  for (int iteration=1; iteration < 1000 && d.fLastResidual() > 1e-25; iteration++) {
  std::cout << "before residual"<<std::endl;
	  residual(g,x);
      d.iterate(g,x);
      std::cout << "iteration "<<iteration<<", Residual norm = "<<std::sqrt(d.fLastResidual())
//                  << ", Distance from solution = "<<std::sqrt((x[0]-1)*(x[0]-1)+(x[1]-1)*(x[1]-1))
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
