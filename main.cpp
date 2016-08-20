#include <iostream>
#include <vector>
#include <cmath>
#include "diiscxx.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace std;

typedef std::vector<double> doubles;

doubles residual(const doubles x) {
  doubles result;
//  result.push_back(2*x[0]-2+400*x[0]*(x[0]*x[0]-x[1])); result.push_back(200*(x[1]-x[0]*x[0])); // Rosenbrock
  result.push_back(x[0]-1);
  result.push_back(x[1]-1);
  return result;
}

doubles hessd(const doubles x) {
  doubles result;
//     result.push_back(2+1200*x[0]*x[0]-400*x[1]); result.push_back(200); // Rosenbrock
  result.push_back(1);
  result.push_back(1);
  return result;
}


using Eigen::MatrixXd;
int main(int argc, char *argv[])
{
  cout << "Test DIIS" << endl;
  doubles amp;amp.push_back(0.99);amp.push_back(0.99);
  doubles res;res.push_back(1);res.push_back(1);
  diis d(amp.size(),res.size(),0,6,1e+6,diis::DIIS);
//  diis DIIS(amp.size(),res.size(),nullptr);
  for (int i=0; i<7& std::fabs(res[0]) > 1e-10 && std::fabs(res[1]) > 1e-10  ; i++) {
    res=residual(amp);
    std::cout <<i<< "before extrapolate amp: "<<amp[0]-1<<", "<<amp[1]-1<<std::endl;
    std::cout <<i<< "before extrapolate res: "<<res[0]<<", "<<res[1]<<std::endl;
    double shift=1;
    std::cout << "amp.data() "<<amp.data()<<std::endl;
    double* a=amp.data();
    d.extrapolate(a,&res[0],nullptr);
    std::cout <<i<< "after  extrapolate amp: "<<amp[0]-1<<", "<<amp[1]-1<<std::endl;
    std::cout <<i<< "after  extrapolate res: "<<res[0]<<", "<<res[1]<<std::endl;
     amp[0]=amp[0]-res[0]/(shift+hessd(amp)[0]);
     amp[1]=amp[1]-res[1]/(shift+hessd(amp)[1]);
    }
  return 0;
}
