#include "Diis.h"
#include <stdexcept>
#include <math.h>

using namespace IterativeSolver;

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
  if (nDim != m_subspaceMatrix.rows()) throw std::logic_error("problem in pruning");
  if (m_Weights.size() != m_subspaceMatrix.rows()) {
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
  for (size_t k=0; k<Coeffs.rows(); k++)
    if (isnan(Coeffs(k))) {
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
static void _Rosenbrock_residual(const ParameterVectorSet & psx, ParameterVectorSet & outputs, std::vector<ParameterScalar> shift=std::vector<ParameterScalar>(), bool append=false) {
  if (not append) outputs.front()->zero();
  (*outputs.front())[0]+=(2*(*psx.front())[0]-2+400*(*psx.front())[0]*((*psx.front())[0]*(*psx.front())[0]-(*psx.front())[1])); (*outputs.front())[1]+=(200*((*psx.front())[1]-(*psx.front())[0]*(*psx.front())[0])); // Rosenbrock
}

static void _Rosenbrock_updater(const ParameterVectorSet & psg, ParameterVectorSet & psc, std::vector<ParameterScalar> shift, bool append=true) {
  if (append) {
      (*psc.front())[0] -=(*psg.front())[0]/700;
      (*psc.front())[1] -=(*psg.front())[1]/200;
    } else
    {
//      xout << "Rosenbrock preconditioner g="<<psg.front()<<std::endl;
      (*psc.front())[0] =-(*psg.front())[0]/700;
      (*psc.front())[1] =-(*psg.front())[1]/200;
//      xout << "Rosenbrock preconditioner c="<<psc.front()<<std::endl;
    }
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
  (*x.front())[0]=(*x.front())[1]=1-difficulty; // initial guess
  xout << "initial guess"<<x[0]<<std::endl;
  xout << "initial guess"<<x<<std::endl;
  bool converged=false;
  for (int iteration=1; iteration < 1000 && not converged; iteration++) {
      _Rosenbrock_residual(x,g);
      optionMap o; //o["weight"]=2;
      converged = d.iterate(g,x,o);
      if (verbosity>0)
        xout << "new x after iterate "<<x.front()<<std::endl;
      if (verbosity>=0)
        xout << "iteration "<<iteration<<", Residual norm = "<<std::sqrt(d.fLastResidual())
             << ", Distance from solution = "<<std::sqrt(((*x.front())[0]-1)*((*x.front())[0]-1)+((*x.front())[1]-1)*((*x.front())[1]-1))
            <<", converged? "<<converged
           <<std::endl;
    }

  if (verbosity>=0) xout << "Test DIIS::solver, difficulty="<<difficulty<<std::endl;
  d.Reset();
  (*x.front())[0]=(*x.front())[1]=1-difficulty; // initial guess
  d.m_verbosity=verbosity;
//  d.solve(g,x);
  xout  << "Distance from solution = "<<std::sqrt(((*x.front())[0]-1)*((*x.front())[0]-1)+((*x.front())[1]-1)*((*x.front())[1]-1)) <<std::endl;

}

#include <random>
#include <chrono>
  struct anharmonic {
    Eigen::MatrixXd m_F;
    double m_gamma;
    size_t m_n;
    anharmonic(){}
    void set(size_t n, double alpha, double gamma)
    {
      m_gamma=gamma;
      m_n=n;
      unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
      std::default_random_engine generator (seed);
      std::uniform_real_distribution<double> distribution(-0.5,0.5);

      m_F.resize(n,n);
      for (size_t j=0; j<n; j++) {
          for (size_t i=0; i<n; i++)
            m_F(i,j)=distribution(generator);
          m_F(j,j) += (j*alpha+0.5);
        }
    }
    SimpleParameterVector guess()
    {
      SimpleParameterVector result(m_n);
      double value=0.3;
      for (size_t k=0; k<m_n; k++) {
        result[k]=value;
        value=-value;
        }
      return result;
    }
  };

  static anharmonic instance;

    static void _anharmonic_residual(const ParameterVectorSet & psx, ParameterVectorSet & outputs, std::vector<ParameterScalar> shift=std::vector<ParameterScalar>(), bool append=false) {
      if (not append) outputs.front()->zero();
      for (size_t i=0; i<instance.m_n; i++) {
          (*outputs.front())[i] = instance.m_gamma*(*psx.front())[i];
          for (size_t j=0; j<instance.m_n; j++)
            (*outputs.front())[i] += instance.m_F(j,i)*(*psx.front())[j];
        }
    }
    static void _anharmonic_preconditioner(const ParameterVectorSet & psg, ParameterVectorSet & psc, std::vector<ParameterScalar> shift=std::vector<ParameterScalar>(), bool append=false) {
      if (append) {
          for (size_t i=0; i<instance.m_n; i++)
            (*psc.front())[i] -= (*psg.front())[i]/instance.m_F(i,i);
        } else {
          for (size_t i=0; i<instance.m_n; i++)
            (*psc.front())[i] =- (*psg.front())[i]/instance.m_F(i,i);
        }
    }
void DIIS::randomTest(size_t sample, size_t n, double alpha, double gamma, DIISmode_type mode)
{

  int nfail=0;
  unsigned int iterations=0, maxIterations=0;
  for (size_t repeat=0; repeat < sample; repeat++) {
      instance.set(n,alpha,gamma);
      DIIS d(&_anharmonic_preconditioner,&_anharmonic_residual);
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
