#include "Diis.h"
#include <stdexcept>

using namespace IterativeSolver;

DIIS::DIIS(ParameterSetTransformation updateFunction, ParameterSetTransformation residualFunction)
  : IterativeSolverBase(updateFunction, residualFunction)
{
  m_orthogonalize = false;
  m_subspaceMatrixResRes = true;
  setOptions();
  Reset();
}

DIIS::~DIIS()
{
}

void DIIS::setOptions(size_t maxDim, double svdThreshold, DIISmode_type mode)
{
  m_maxDim = maxDim;
  m_svdThreshold = svdThreshold;
  m_DIISmode = mode;

  Reset();
  if (m_verbosity>1)
    xout << "m_DIISmode set to "<<m_DIISmode<<" in setOptions"<<std::endl;
  if (m_verbosity>1)
    xout << "m_maxDim set to "<<m_maxDim<<" in setOptions"<<std::endl;
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
  m_LastResidualNormSq = m_subspaceMatrix(nDim-1,nDim-1);

  m_Weights.push_back(weight);

  int worst=0, best=0;
  for (size_t i=0; i<nDim; i++) {
      if (m_subspaceMatrix(i,i) > m_subspaceMatrix(worst,worst)) worst=i;
      if (m_subspaceMatrix(i,i) < m_subspaceMatrix(best,best)) best=i;
  }
  double fBaseScale = std::sqrt(m_subspaceMatrix(worst,worst)*m_subspaceMatrix(best,best));

  if (nDim > m_maxDim) { // prune away the worst/oldest vector. Algorithm to be done properly yet
      size_t prune = worst;
      if (true || prune == nDim-1) { // the current vector is the worst, so delete the oldest
//          xout << "m_dataOfBirth: "; for (auto b=m_dateOfBirth.begin(); b!=m_dateOfBirth.end(); b++) xout <<(*b); xout<<std::endl;
          prune = std::min_element(m_dateOfBirth.begin(),m_dateOfBirth.end())-m_dateOfBirth.begin();
      }
//      xout << "prune="<<prune<<std::endl;
      deleteVector(prune);
      m_Weights.erase(m_Weights.begin()+prune);
      nDim--;
  }
  if (nDim != m_subspaceMatrix.rows()) throw std::logic_error("problem in pruning");
  if (m_Weights.size() != m_subspaceMatrix.rows()) throw std::logic_error("problem in pruning weights");



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
  Rhs[nDim] = -weight;
//  xout << "B:"<<std::endl<<B<<std::endl;
//  xout << "Rhs:"<<std::endl<<Rhs<<std::endl;

  // invert the system, determine extrapolation coefficients.
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(B, Eigen::ComputeThinU | Eigen::ComputeThinV);
  svd.setThreshold(1e-10);
  Coeffs=svd.solve(Rhs);
  if (m_verbosity>1) xout << "Combination of iteration vectors: "<<Coeffs.transpose()<<std::endl;

  m_LastAmplitudeCoeff = Coeffs[nDim-1];
  //  xout << "m_solutions.front(): "<<m_solutions.front()<<std::endl;
  for (size_t l=0; l<residual.size(); l++) residual[l]->zero();
  //  xout << "m_solutions.front(): "<<m_solutions.front()<<std::endl;
  for (size_t l=0; l<solution.size(); l++) solution[l]->zero();
  //  xout << "residual at "<<&residual[0]<<std::endl;
  //  xout << "solution at "<<&solution[0]<<std::endl;
  //  xout << "m_solutions.front() at "<<&m_solutions.front()[0]<<std::endl;
  //  xout << "m_solutions.size(): "<<m_solutions.size()<<std::endl;
  //  xout << "m_solutions.front(): "<<m_solutions.front()<<std::endl;
  //  xout << "m_solutions.front()[0]: "<<m_solutions.front()[0]<<std::endl;
  size_t k=0;
  for (size_t l=0; l<other.size(); l++) other[l]->zero();
  for (size_t kk=0; kk<m_residuals.size(); kk++) {
      if (m_residuals[kk].m_active.front()){
//            xout << "extrapolation k="<<k<<std::endl;
//            xout << "add solution "<<kk<<m_solutions[kk].front()<<" with factor "<<Coeffs[k]<<std::endl;
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
  if (not append) psc.front()->zero();
  (*psc.front())[0] -=(*psg.front())[0]/700;
  (*psc.front())[1] -=(*psg.front())[1]/200;
}

void DIIS::test(int verbosity,
                size_t maxDim, double acceptanceThreshold, DIISmode_type mode, double difficulty)
{
  SimpleParameterVector xx(2);
  SimpleParameterVector gg(2);
  ParameterVectorSet x; x.push_back(&xx);
  ParameterVectorSet g; g.push_back(&gg);
  DIIS d(&_Rosenbrock_updater,&_Rosenbrock_residual);
  d.setOptions(maxDim, acceptanceThreshold, mode);

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

