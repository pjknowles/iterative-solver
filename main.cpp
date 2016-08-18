#include <iostream>
#include <vector>
#include <cmath>
#include "diiscxx.h"
#include <Eigen/Dense>

using namespace std;

typedef std::vector<double> doubles;

doubles residual(const doubles x) {
  doubles result;
  result.push_back(2*x[0]-2+400*x[0]*(x[0]*x[0]-x[1]));
  result.push_back(200*(x[1]-x[0]*x[0]));
  return result;
}

using Eigen::MatrixXd;
int main(int argc, char *argv[])
{
  cout << "Hello World!" << endl;
  MatrixXd m(2,2);
  m(0,0) = 3;
  m(1,0) = 2.5;
  m(0,1) = -1;
  m(1,1) = m(1,0) + m(0,1);
  std::cout << m << std::endl;
  doubles amp(2,0.99);
  diiscxx d(amp.size());
  doubles resid({1,1});
  for (int i=0; i<10&& std::fabs(resid[0]) > 1e-10 && std::fabs(resid[1]) > 1e-10  ; i++) {
    resid=residual(amp);
    std::cout <<i<< "amp: "<<amp[0]-1<<", "<<amp[1]-1<<std::endl;
    std::cout <<i<< "resid: "<<resid[0]<<", "<<resid[1]<<std::endl;
    double shift=1;
     amp[0]=amp[0]-resid[0]/(shift+2+1200*amp[0]*amp[0]-400*amp[1]);
     amp[1]=amp[1]-resid[1]/(shift+200);
    }
  return 0;
}
