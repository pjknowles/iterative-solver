#include <iostream>
#include <vector>
#include <cmath>
#include "Diis.h"

using namespace std;
using namespace IterativeSolver;

typedef std::vector<double> doubles;

bool rosenbrock=true;

void globalsum(double* buffer, size_t length)
{
  std::cout << "globalsum"; for (size_t k=0; k<length; k++) std::cout << " "<<buffer[k]; std:cout << endl;
}

doubles residual(const doubles x) {
  doubles result;
  if (rosenbrock) {
  result.push_back(2*x[0]-2+400*x[0]*(x[0]*x[0]-x[1])); result.push_back(200*(x[1]-x[0]*x[0])); // Rosenbrock
    } else {
  result.push_back(x[0]-1);
  result.push_back(x[1]-1);
    }
  return result;
}

doubles hessd(const doubles x) {
  doubles result;
  if (rosenbrock) {
     result.push_back(2+1200*x[0]*x[0]-400*x[1]); result.push_back(200); // Rosenbrock
    } else {
  result.push_back(1);
  result.push_back(1);
    }
  return result;
}

void simple() {
  std::vector<double> x(2);
  std::vector<double> g(2);
  x[0]=x[1]=0.9; // initial guess
  std::vector<double> diag(2); diag[0]=700; diag[1]=200; // preconditioner
  std::vector<size_t> lengths; lengths.push_back(g.size()); lengths.push_back(x.size());
  Diis d(lengths);
//  d.setGlobalSum(&globalsum);
//  auto gs = [](double* buffer, size_t length){std::cout << "hello"<<std::endl;};
//  d.setGlobalSum(&gs);
  d.addPreconditioner(&diag[0],0,true);
  for (int iteration=1; iteration < 1000 && d.fLastResidual() > 1e-25; iteration++) {
      g[0]=2*x[0]-2+400*x[0]*(x[0]*x[0]-x[1]); g[1]=200*(x[1]-x[0]*x[0]); // Rosenbrock function gradient
      d.iterate(&g[0],&x[0]);
      std::cout << "iteration "<<iteration<<", Residual norm = "<<std::sqrt(d.fLastResidual())
                  << ", Distance from solution = "<<std::sqrt((x[0]-1)*(x[0]-1)+(x[1]-1)*(x[1]-1))<<std::endl;
    }
}

int main(int argc, char *argv[])
{
  cout << "Test DIIS" << endl;
  doubles amp;amp.push_back(0.9);amp.push_back(0.9);
  doubles res;res.push_back(1);res.push_back(1);
  std::vector<size_t> lengths; lengths.push_back(res.size()); lengths.push_back(amp.size());
  Diis d(lengths,6,1e+6,Diis::DIIS);
  bool iterate=true;
  int verbosity=0;
  d.setVerbosity(verbosity);
  double shift=1;
  d.addPreconditioner(&hessd(amp)[0],shift,true);
  std::cout << "jacobian diagonals "<<hessd(amp)[0]<<" "<<hessd(amp)[1]<<std::endl;
//  d.setBufferSize(2);
  for (int i=1; i<10000 && d.fLastResidual() > 1e-25; i++) {
    res=residual(amp);
    if (verbosity>0) {
    std::cout <<i<< "before extrapolate amp: "<<amp[0]-1<<", "<<amp[1]-1<<std::endl;
    std::cout <<i<< "before extrapolate res: "<<res[0]<<", "<<res[1]<<std::endl;
      }
    if (iterate) {
        d.iterate(&res[0],&amp[0]);
        std::cout << "Iteration "<<i<<", Residual norm = "<<std::sqrt(d.fLastResidual())
                  << ", Distance from solution = "<<std::sqrt((amp[0]-1)*(amp[0]-1)+(amp[1]-1)*(amp[1]-1))<<std::endl;
      } else {
//    d.extrapolate({&amp[0],&res[0]});
    std::vector<double*> vectors;vectors.push_back(&res[0]);vectors.push_back(&amp[0]);
    d.extrapolate(vectors);
    std::cout << "Iteration "<<i<<", Residual norm = "<<std::sqrt(d.fLastResidual())
              << ", Distance from solution = "<<std::sqrt((amp[0]-1)*(amp[0]-1)+(amp[1]-1)*(amp[1]-1))<<std::endl;
    if (verbosity>0) {
    std::cout <<i<< "fLastResidual="<<d.fLastResidual()<<std::endl;
    std::cout <<i<< "fLastCoeff="<<d.fLastCoeff()<<std::endl;
    std::cout <<i<< "after  extrapolate amp: "<<amp[0]-1<<", "<<amp[1]-1<<std::endl;
    std::cout <<i<< "after  extrapolate res: "<<res[0]<<", "<<res[1]<<std::endl;
      }
    if (false) {
     amp[0]=amp[0]-res[0]/(shift+hessd(amp)[0]);
     amp[1]=amp[1]-res[1]/(shift+hessd(amp)[1]);
      } else {
        d.update(&res[0],&amp[0]);
      }

      }
    }
  std::cout << "try simple version"<<std::endl;
  simple();
  return 0;
}
