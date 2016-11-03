#include "ISDiis.h"
#include <stdexcept>
#include <math.h>

using namespace LinearAlgebra;

DIIS::DIIS(const ParameterSetTransformation residualFunction, const ParameterSetTransformation updateFunction)
  : IterativeSolverBase(residualFunction, updateFunction)
  , m_maxDim(6)
  , m_svdThreshold(1e-10)
{
  m_orthogonalize = false;
  setMode(DIISmode);
  Reset();
}

DIIS::~DIIS()
{
}

void DIIS::setMode(DIISmode_type mode)
{
  m_DIISmode = mode;
  m_subspaceMatrixResRes = mode!=KAINmode;
  m_preconditionResiduals = mode==KAINmode;

  Reset();
  if (m_verbosity>1)
    xout << "m_DIISmode set to "<<m_DIISmode<<std::endl;
}

void DIIS::Reset()
{
  m_LastResidualNormSq=1e99; // so it can be tested even before extrapolation is done
  m_lastVectorIndex=0;
  while (m_subspaceMatrix.rows()>0) deleteVector(0); ;
  m_residuals.clear();
  m_solutions.clear();
  m_others.clear();
  m_Weights.clear();
}


void DIIS::extrapolate(ParameterVectorSet & residual, ParameterVectorSet & solution, ParameterVectorSet & other, const optionMap options)
{
  //	  xout << "Enter DIIS::extrapolate"<<std::endl;
  //	  xout << "residual : "<<residual<<std::endl;
  //	  xout << "solution : "<<solution<<std::endl;
  m_updateShift.clear();m_updateShift.push_back(-(1+std::numeric_limits<double>::epsilon())*m_subspaceMatrix(0,0));
  double weight=options.count("weight") ? (options.find("weight")->second) : 1.0;
  if (m_maxDim <= 1 || m_DIISmode == disabled) return;

  if (residual.size() > 1) throw std::logic_error("DIIS does not handle multiple solutions");
  m_lastVectorIndex=m_residuals.size()-1;

//  if (m_subspaceMatrix.rows() < 9) {
//      xout << "m_subspaceMatrix on entry to DIIS::extrapolate"<<std::endl<<m_subspaceMatrix<<std::endl;
//  }
  size_t nDim = m_subspaceMatrix.rows();
  m_LastResidualNormSq = std::fabs(m_subspaceMatrix(nDim-1,nDim-1));
//  xout << "m_LastResidualNormSq "<<m_LastResidualNormSq<<std::endl;

  m_Weights.push_back(weight);

  Eigen::ArrayXd d = m_subspaceMatrix.diagonal().array().abs();
  int worst=0, best=0;
  for (size_t i=0; i<nDim; i++) {
      if (d(i) > d(worst)) worst=i;
      if (d(i) < d(best)) best=i;
  }
  double fBaseScale = std::sqrt(d(worst)*d(best));

  while (nDim > m_maxDim) { // prune away the worst/oldest vector. Algorithm to be done properly yet
      size_t prune = worst;
      if (true || prune == nDim-1) { // the current vector is the worst, so delete the oldest
//          xout << "m_dateOfBirth: "; for (auto b=m_dateOfBirth.begin(); b!=m_dateOfBirth.end(); b++) xout <<(*b); xout<<std::endl;
          prune = std::min_element(m_dateOfBirth.begin(),m_dateOfBirth.end())-m_dateOfBirth.begin();
      }
//      xout << "prune="<<prune<<std::endl;
//  xout << "nDim="<<nDim<<", m_Weights.size()="<<m_Weights.size()<<std::endl;
      deleteVector(prune);
      m_Weights.erase(m_Weights.begin()+prune);
      nDim--;
//  xout << "nDim="<<nDim<<", m_Weights.size()="<<m_Weights.size()<<std::endl;
  }
  if (nDim != (size_t)m_subspaceMatrix.rows()) throw std::logic_error("problem in pruning");
  if (m_Weights.size() != (size_t)m_subspaceMatrix.rows()) {
      xout << "nDim="<<nDim<<", m_Weights.size()="<<m_Weights.size()<<std::endl;
      throw std::logic_error("problem after pruning weights");
    }



//  if (m_subspaceMatrix.rows() < 9) {
//      xout << "m_subspaceMatrix"<<std::endl<<m_subspaceMatrix<<std::endl;
//  }

  // build actual DIIS system for the subspace used.
  Eigen::VectorXd
      Rhs(nDim+1),
      Coeffs(nDim+1);
  Eigen::MatrixXd
      B(nDim+1, nDim+1);

  // Factor out common size scales from the residual dots.
  // This is done to increase numerical stability for the case when _all_
  // residuals are very small.
      B.block(0,0,nDim,nDim) = m_subspaceMatrix/fBaseScale;
      Rhs.head(nDim) = Eigen::VectorXd::Zero(nDim);

  // make Lagrange/constraint lines.
  for ( size_t i = 0; i < nDim; ++ i )
      B(i, nDim) =  B(nDim, i) = -m_Weights[i];
  B(nDim, nDim) = 0.0;
  Rhs[nDim] = -1;
//  xout << "B:"<<std::endl<<B<<std::endl;
//  xout << "Rhs:"<<std::endl<<Rhs<<std::endl;

  // invert the system, determine extrapolation coefficients.
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(B, Eigen::ComputeThinU | Eigen::ComputeThinV);
  svd.setThreshold(m_svdThreshold);
  Coeffs=svd.solve(Rhs).head(nDim);
  m_LastAmplitudeCoeff = Coeffs[nDim-1];
  if (m_verbosity>1) xout << "Combination of iteration vectors: "<<Coeffs.transpose()<<std::endl;
  for (size_t k=0; k<(size_t)Coeffs.rows(); k++)
    if (std::isnan(Coeffs(k))) {
        xout << "B:"<<std::endl<<B<<std::endl;
        xout << "Rhs:"<<std::endl<<Rhs<<std::endl;
        xout << "Combination of iteration vectors: "<<Coeffs.transpose()<<std::endl;
        throw std::overflow_error("NaN detected in DIIS submatrix solution");
      }
  residual.zero();
  solution.zero();
  size_t k=0;
  for (size_t l=0; l<other.size(); l++) other[l]->zero();
  for (size_t kk=0; kk<m_residuals.size(); kk++) {
      if (m_residuals[kk].m_active.front()){
      residual.front()->axpy(Coeffs[k],m_residuals[kk].front());
      solution.front()->axpy(Coeffs[k],m_solutions[kk].front());
      for (size_t l=0; l<other.size(); l++)
        other[l]->axpy(Coeffs[k],m_others[kk][l]);
      k++;
      }
    }
  if (m_verbosity>2) {
      xout << "DIIS.extrapolate() final extrapolated solution: "<<solution.front()<<std::endl;
      xout << "DIIS.extrapolate() final extrapolated residual: "<<residual.front()<<std::endl;
    }
}

// testing code below here
#include "SimpleParameterVector.h"
static void _Rosenbrock_residual(const ParameterVectorSet & psx, ParameterVectorSet & outputs, std::vector<scalar> shift=std::vector<scalar>(), bool append=false) {
  size_t n=2;
  std::vector<scalar> psxk(n);
  std::vector<scalar> output(n);

  if (not append) outputs.front()->zero();
  outputs.front()->get(&(psxk[0]),n,0);
  psx.front()->get(&(psxk[0]),n,0);
  output[0]+=(2*psxk[0]-2+400*psxk[0]*(psxk[0]*psxk[0]-psxk[1])); output[1]+=(200*(psxk[1]-psxk[0]*psxk[0])); // Rosenbrock
  outputs.front()->put(&(psxk[0]),n,0);
}

static void _Rosenbrock_updater(const ParameterVectorSet & psg, ParameterVectorSet & psc, std::vector<scalar> shift, bool append=true) {
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

void DIIS::test(int verbosity,
                size_t maxDim, double svdThreshold, DIISmode_type mode, double difficulty)
{
  SimpleParameterVector xx(2);
  SimpleParameterVector gg(2);
  ParameterVectorSet x; x.push_back(&xx);
  ParameterVectorSet g; g.push_back(&gg);
  DIIS d(&_Rosenbrock_residual,&_Rosenbrock_updater);
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
      converged = d.iterate(g,x,o);
      x.front()->get(&xxx[0],2,0);
      if (verbosity>2)
        xout << "new x after iterate "<<x.front()<<std::endl;
      if (verbosity>=0)
        xout << "iteration "<<iteration<<", Residual norm = "<<std::sqrt(d.fLastResidual())
             << ", Distance from solution = "<<std::sqrt((xxx[0]-1)*(xxx[0]-1)+(xxx[1]-1)*(xxx[1]-1))
            <<", converged? "<<converged
           <<std::endl;
    }

  if (verbosity>=0) xout << "Test DIIS::solver, difficulty="<<difficulty<<std::endl;
  d.Reset();
  xxx[0]=xxx[1]=1-difficulty; // initial guess
  xx.put(&xxx[0],2,0);
  d.m_verbosity=verbosity;
//  d.solve(g,x);
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

  static anharmonic instance;

    static void _anharmonic_residual(const ParameterVectorSet & psx, ParameterVectorSet & outputs, std::vector<scalar> shift=std::vector<scalar>(), bool append=false) {
      std::vector<scalar> psxk(instance.m_n);
      std::vector<scalar> output(instance.m_n);
      psx.front()->get(&(psxk[0]),instance.m_n,0);
      if (append)
        outputs.front()->get(&(output[0]),instance.m_n,0);
      else
        outputs.front()->zero();

      for (size_t i=0; i<instance.m_n; i++) {
          output[i] = instance.m_gamma*psxk[i];
          for (size_t j=0; j<instance.m_n; j++)
            output[i] += instance.m_F(j,i)*psxk[j];
        }
      outputs.front()->put(&output[0],instance.m_n,0);
    }
    static void _anharmonic_preconditioner(const ParameterVectorSet & psg, ParameterVectorSet & psc, std::vector<scalar> shift=std::vector<scalar>(), bool append=false) {
      std::vector<scalar> psck(instance.m_n);
      std::vector<scalar> psgk(instance.m_n);
      psg.front()->get(&psgk[0],instance.m_n,0);
      if (append) {
          psc.front()->get(&psck[0],instance.m_n,0);
          for (size_t i=0; i<instance.m_n; i++)
            psck[i] -= psgk[i]/instance.m_F(i,i);
        } else {
          for (size_t i=0; i<instance.m_n; i++)
            psck[i] =- psgk[i]/instance.m_F(i,i);
        }
      psc.front()->put(&psck[0],instance.m_n,0);
    }
void DIIS::randomTest(size_t sample, size_t n, double alpha, double gamma, DIISmode_type mode)
{

  int nfail=0;
  unsigned int iterations=0, maxIterations=0;
  for (size_t repeat=0; repeat < sample; repeat++) {
      instance.set(n,alpha,gamma);
      DIIS d(&_anharmonic_residual,&_anharmonic_preconditioner);
      d.setMode(mode);
      d.m_verbosity=-1;
      d.m_maxIterations=100000;
      SimpleParameterVector gg(n); ParameterVectorSet g; g.push_back(&gg);
      SimpleParameterVector xx=instance.guess(); ParameterVectorSet x; x.push_back(&xx);
      if (not d.solve(g,x)) nfail++;
      iterations+=d.iterations();
      if (maxIterations<d.iterations())
        maxIterations=d.iterations();
    }
  xout << "sample="<<sample<<", n="<<n<<", alpha="<<alpha<<", gamma="<<gamma<<", average iterations="<<iterations/sample<<", maximum iterations="<<maxIterations<<", nfail="<<nfail<<std::endl;
}
