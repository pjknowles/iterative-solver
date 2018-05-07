#include "ISDiis.h"
#include <stdexcept>
#include <memory>
#include <math.h>
#include <Eigen/Jacobi>
#include <cmath>

using namespace LinearAlgebra;


// testing code below here
#include "SimpleParameterVector.h"
using scalar=double;
template<>
void DIIS<scalar>::test(int verbosity,
                size_t maxDim, double svdThreshold, DIISmode_type mode, double difficulty)
{
static struct {
  void operator()(const vectorSet<scalar> & psx, vectorSet<scalar> & outputs) const {
    size_t n=2;
    std::vector<scalar> psxk(n);
    std::vector<scalar> output(n);

    psx.front()->get(&(psxk[0]),n,0);
    output[0]=(2*psxk[0]-2+400*psxk[0]*(psxk[0]*psxk[0]-psxk[1])); output[1]=(200*(psxk[1]-psxk[0]*psxk[0])); // Rosenbrock
    outputs.front()->put(&(output[0]),n,0);
  }
} _Rosenbrock_residual;

static struct {
  void operator()(const vectorSet<scalar> & psg, vectorSet<scalar> & psc, std::vector<scalar> shift, bool append=true) const {
    size_t n=2;
    std::vector<scalar> psck(n);
    std::vector<scalar> psgk(n);
    psg.front()->get(&psgk[0],n,0);
    if (append) {
        psc.front()->get(&psck[0],n,0);
        psck[0] -=psgk[0]/700;
        psck[1] -=psgk[1]/200;
      } else
      {
        psck[0] =-psgk[0]/700;
        psck[1] =-psgk[1]/200;
      }
    psc.front()->put(&psck[0],n,0);
  }
} _Rosenbrock_updater;
  SimpleParameterVector xx(2);
  SimpleParameterVector gg(2);
  vectorSet<scalar> x; x.push_back(std::shared_ptr<SimpleParameterVector>(&xx));
  vectorSet<scalar> g; g.push_back(std::shared_ptr<SimpleParameterVector>(&gg));
  DIIS d;
  d.m_maxDim=maxDim;
  d.m_svdThreshold=svdThreshold;
  d.setMode(mode);

  if (verbosity>=0) xout << "Test DIIS::iterate, difficulty="<<difficulty<<std::endl;
  d.Reset();
  d.m_verbosity=verbosity-1;
  std::vector<scalar> xxx(2);
  xxx[0]=xxx[1]=1-difficulty; // initial guess
  xx.put(&xxx[0],2,0);
  xout << "initial guess"<<x[0]<<std::endl;
  xout << "initial guess"<<x<<std::endl;
  bool converged=false;
  for (int iteration=1; iteration < 1000 && not converged; iteration++) {
      _Rosenbrock_residual(x,g);
      optionMap o; //o["weight"]=2;
      d.interpolate(x,g,o);
      std::vector<scalar> shift; shift.push_back(1e-10);
      _Rosenbrock_updater(x,g,shift);
      converged = d.finalize(x,g);
      x.front()->get(&xxx[0],2,0);
      if (verbosity>2)
        xout << "new x after iterate "<<x.front()<<std::endl;
      if (verbosity>=0)
        xout << "iteration "<<iteration<<", Residual norm = "<<std::sqrt(d.fLastResidual())
             << ", Distance from solution = "<<std::sqrt((xxx[0]-1)*(xxx[0]-1)+(xxx[1]-1)*(xxx[1]-1))
            <<", converged? "<<converged
           <<std::endl;
    }

  x.front()->get(&xxx[0],2,0);
  xout   << "Distance from solution = "<<std::sqrt((xxx[0]-1)*(xxx[0]-1)+(xxx[1]-1)*(xxx[1]-1));

}

#include <cstdlib>
struct anharmonic {
  Eigen::MatrixXd m_F;
  double m_gamma;
  size_t m_n;
  anharmonic(){}
  void set(size_t n, double alpha, double gamma)
  {
    m_gamma=gamma;
    m_n=n;

    m_F.resize(n,n);
    for (size_t j=0; j<n; j++) {
        for (size_t i=0; i<n; i++)
          m_F(i,j)=-0.5 + (((double)rand())/RAND_MAX);
        m_F(j,j) += (j*alpha+0.5);
      }
  }
  SimpleParameterVector guess()
  {
    std::vector<scalar> r(m_n);
    SimpleParameterVector result(m_n);
    double value=0.3;
    for (size_t k=0; k<m_n; k++) {
        r[k]=value;
        value=-value;
      }
    result.put(&r[0],m_n,0);
    return result;
  }
};

//static anharmonic instance;

//static struct : IterativeSolverBase::ParameterSetTransformation {
//  void operator()(const vectorSet<scalar> & psx, vectorSet<scalar> & outputs, std::vector<scalar> shift=std::vector<scalar>(), bool append=false) const override {
//    std::vector<scalar> psxk(instance.m_n);
//    std::vector<scalar> output(instance.m_n);
//    psx.front()->get(&(psxk[0]),instance.m_n,0);
//    if (append)
//      outputs.front()->get(&(output[0]),instance.m_n,0);
//    else
//      outputs.front()->zero();

//    for (size_t i=0; i<instance.m_n; i++) {
//        output[i] = instance.m_gamma*psxk[i];
//        for (size_t j=0; j<instance.m_n; j++)
//          output[i] += instance.m_F(j,i)*psxk[j];
//      }
//    outputs.front()->put(&output[0],instance.m_n,0);
//  }
//} _anharmonic_residual;
//static struct : IterativeSolverBase::ParameterSetTransformation {
//  void operator()(const vectorSet<scalar> & psg, vectorSet<scalar> & psc, std::vector<scalar> shift=std::vector<scalar>(), bool append=false) const override {
//    std::vector<scalar> psck(instance.m_n);
//    std::vector<scalar> psgk(instance.m_n);
//    psg.front()->get(&psgk[0],instance.m_n,0);
//    if (append) {
//        psc.front()->get(&psck[0],instance.m_n,0);
//        for (size_t i=0; i<instance.m_n; i++)
//          psck[i] -= psgk[i]/instance.m_F(i,i);
//      } else {
//        for (size_t i=0; i<instance.m_n; i++)
//          psck[i] =- psgk[i]/instance.m_F(i,i);
//      }
//    psc.front()->put(&psck[0],instance.m_n,0);
//  }
//} _anharmonic_preconditioner;
//void DIIS::randomTest(size_t sample, size_t n, double alpha, double gamma, DIISmode_type mode)
//{

//  int nfail=0;
//  unsigned int iterations=0, maxIterations=0;
//  for (size_t repeat=0; repeat < sample; repeat++) {
//      instance.set(n,alpha,gamma);
//      DIIS d(_anharmonic_residual,_anharmonic_preconditioner);
//      d.setMode(mode);
//      d.m_verbosity=-1;
//      d.m_maxIterations=100000;
//      SimpleParameterVector gg(n); vectorSet<scalar> g; g.push_back(std::shared_ptr<SimpleParameterVector>(&gg));
//      SimpleParameterVector xx=instance.guess(); vectorSet<scalar> x; x.push_back(std::shared_ptr<SimpleParameterVector>(&xx));
//      if (not d.solve(g,x)) nfail++;
//      iterations+=d.iterations();
//      if (maxIterations<d.iterations())
//        maxIterations=d.iterations();
//    }
//  xout << "sample="<<sample<<", n="<<n<<", alpha="<<alpha<<", gamma="<<gamma<<", average iterations="<<iterations/sample<<", maximum iterations="<<maxIterations<<", nfail="<<nfail<<std::endl;
//}
