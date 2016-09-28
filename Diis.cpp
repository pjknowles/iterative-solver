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

void DIIS::setOptions(size_t maxDim, double acceptanceThreshold, DIISmode_type mode)
{
  m_maxDim = maxDim;
  m_acceptanceThreshold = acceptanceThreshold;
  m_DIISmode = mode;

  Reset();
  if (m_verbosity>1)
    xout << "m_DIISmode set to "<<m_DIISmode<<" in setOptions"<<std::endl;
  if (m_verbosity>1)
    xout << "m_maxDim set to "<<m_maxDim<<" in setOptions"<<std::endl;
}

void DIIS::Reset()
{
  m_iNext=0;
  m_iVectorAge.resize(0);
  m_LastResidualNormSq=1e99; // so it can be tested even before extrapolation is done
}

void DIIS::LinearSolveSymSvd(Eigen::VectorXd& Out, const Eigen::MatrixXd& Mat, const Eigen::VectorXd& In, unsigned int nDim, double Thr)
{
  Eigen::VectorXd Ews(nDim);
  Eigen::VectorXd Xv(nDim); // input vectors in EV basis.
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
  es.compute(-Mat);
  Ews = es.eigenvalues();

  if (m_verbosity > 1) {
      xout << "diis::LinearSolveSymSvd"<<std::endl;
      xout << "Mat="<<Mat<<std::endl;
      xout << "Ews="<<Ews<<std::endl;
      xout << "In="<<In<<std::endl;
      xout << "es.eigenvectors()="<<es.eigenvectors()<<std::endl;
    }

  Xv = In.transpose()*es.eigenvectors();
  //    xout << "Xv="<<Xv<<std::endl;
  for (size_t iEw = 0; iEw != nDim; ++ iEw)
    if (std::abs(Ews(iEw)) >= Thr)
      Xv(iEw) /= -Ews(iEw);
  // ^- note that this only screens by absolute value.
  // no positive semi-definiteness is assumed!
    else
      Xv(iEw) = 0.;
  if (m_verbosity > 1) xout << "Xv="<<Xv<<std::endl;
  Out = es.eigenvectors() * Xv;
  Out.conservativeResize(Out.size()-1,1);
  if (m_verbosity > 1) xout << "Out="<<Out<<std::endl;
}



void DIIS::extrapolate(ParameterVectorSet & residual, ParameterVectorSet & solution, ParameterVectorSet & other, std::string options)
{
  //	  xout << "Enter DIIS::extrapolate"<<std::endl;
  //	  xout << "residual : "<<residual<<std::endl;
  //	  xout << "solution : "<<solution<<std::endl;
  double weight=1.0;
  size_t pos=options.find("weight=");
  if (pos != std::string::npos)  throw std::logic_error("parsing of weight not implemented yet"); //FIXME
  if (m_maxDim <= 1 || m_DIISmode == disabled) return;

  if (residual.size() > 1) throw std::logic_error("DIIS does not handle multiple solutions");

  double fThisResidualDot = residual.front()->dot(residual.front());
  m_LastResidualNormSq = fThisResidualDot;

  if (m_verbosity>1)
    xout << "m_iNext="<<m_iNext<<", fThisResidualDot="<<fThisResidualDot<<", m_acceptanceThreshold="<<m_acceptanceThreshold<<std::endl;
  if ( m_iNext == 0 && fThisResidualDot > m_acceptanceThreshold )
    // current vector is to be considered too wrong to be useful for DIIS
    // purposes. Don't store it.
    return;

  uint iThis = m_iNext;
  if (m_verbosity > 1) xout<< "iThis=m_iNext "<<m_iNext<<std::endl;
  assert(iThis < m_maxDim);
  if (iThis >= m_iVectorAge.size()) m_iVectorAge.resize(iThis+1);
  if (iThis >= m_ErrorMatrix.cols()) m_ErrorMatrix.conservativeResize(iThis+1,iThis+1);
  m_ErrorMatrix(iThis,iThis)=fThisResidualDot;
  for ( uint i = 0; i < m_iVectorAge.size(); ++ i )
    m_iVectorAge[i] += 1;
  m_iVectorAge[iThis] = 0;
  if (m_verbosity>1) {
      xout << "iVectorAge:";
      for (std::vector<uint>::const_iterator a=m_iVectorAge.begin(); a!=m_iVectorAge.end(); a++)
        xout << " "<<*a;
      xout << std::endl;
    }

  // find set of vectors actually used in the current run and
  // find their common size scale.
  uint
      nDim;
  std::vector<uint> iUsedVecs( m_iVectorAge.size() + 1 );
  // ^- note: this is the DIIS dimension--the actual matrices and vectors have
  // dimension nDim+1 due to the Lagrange-Multipliers! ?? PJK ??
  double
      fBaseScale;
  FindUsefulVectors(&iUsedVecs[0], nDim, fBaseScale, iThis);
  if (m_verbosity>1) {
      xout << "iUsedVecs:";
      for (std::vector<uint>::const_iterator a=iUsedVecs.begin(); a!=iUsedVecs.end(); a++)
        xout << " "<<*a;
      xout << std::endl;
    }
  // transform iThis into a relative index.
  for ( uint i = 0; i < nDim; ++ i )
    if ( iThis == iUsedVecs[i] ) {
        iThis = i;
        break;
      }
  // store current vectors into their position in the array mapped by iUsedVecs
  if (iUsedVecs[iThis] != m_residuals.size()-1)
    {
      m_solutions[iUsedVecs[iThis]]=solution;m_solutions.pop_back();
      m_residuals[iUsedVecs[iThis]]=residual;m_residuals.pop_back();
      m_others[iUsedVecs[iThis]]=other;m_others.pop_back();
      m_lastVectorIndex=iUsedVecs[iThis];
    }
  //  xout << "after replacing expansion vector"<<std::endl;
  //  for (auto vs=m_solutions.begin(); vs!=m_solutions.end(); vs++)
  //      xout << "solution vector"<<*vs<<std::endl;

  if (m_Weights.size()<=iUsedVecs[iThis]) m_Weights.resize(iUsedVecs[iThis]+1);
  m_Weights[iUsedVecs[iThis]] = weight;

  // go through previous residual vectors and form the dot products with them
  std::vector<double> ResDot;
  for (std::vector<ParameterVectorSet>::iterator p=m_residuals.begin(); p!=m_residuals.end(); p++)
    ResDot.push_back(residual.front()->dot(p->front()));
  if (iThis >= nDim) ResDot.resize(iThis+1);
  ResDot[iThis] = fThisResidualDot;


  // update resident error matrix with new residual-dots
  for ( uint i = 0; i < nDim; ++ i ) {
      m_ErrorMatrix(iUsedVecs[i], iUsedVecs[iThis]) = ResDot[i];
      m_ErrorMatrix(iUsedVecs[iThis], iUsedVecs[i]) = ResDot[i];
    }

  if (m_subspaceMatrix.rows() < 9) {
      xout << "m_ErrorMatrix"<<std::endl<<m_ErrorMatrix<<std::endl;
      xout << "m_subspaceMatrix"<<std::endl<<m_subspaceMatrix<<std::endl;
  }

  // build actual DIIS system for the subspace used.
  Eigen::VectorXd
      Rhs(nDim+1),
      Coeffs(nDim+1);
  Eigen::MatrixXd
      B(nDim+1, nDim+1);

  // Factor out common size scales from the residual dots.
  // This is done to increase numerical stability for the case when _all_
  // residuals are very small.
  for ( uint nRow = 0; nRow < nDim; ++ nRow )
    for ( uint nCol = 0; nCol < nDim; ++ nCol )
      B(nRow, nCol) = m_ErrorMatrix(iUsedVecs[nRow], iUsedVecs[nCol])/fBaseScale;

  // make Lagrange/constraint lines.
  for ( uint i = 0; i < nDim; ++ i ) {
      B(i, nDim) = -m_Weights[iUsedVecs[i]];
      B(nDim, i) = -m_Weights[iUsedVecs[i]];
      Rhs[i] = 0.0;
    }
  B(nDim, nDim) = 0.0;
  Rhs[nDim] = -weight;

  // invert the system, determine extrapolation coefficients.
  LinearSolveSymSvd(Coeffs, B, Rhs, nDim+1, 1.0e-10);
  if (m_verbosity>1) xout << "Combination of iteration vectors: "<<Coeffs.transpose()<<std::endl;

  // Find a storage place for the vector in the next round. Either
  // an empty slot or the oldest vector.
  uint iOldestAge = m_iVectorAge[0];
  m_iNext = 0;
  if (m_iVectorAge.size() < m_maxDim) {
      m_iNext = m_iVectorAge.size(); m_iVectorAge.push_back(0);
    } else
    for ( uint i = m_iVectorAge.size(); i != 0; -- i ){
        if ( iOldestAge <= m_iVectorAge[i-1] ) {
            iOldestAge = m_iVectorAge[i-1];
            m_iNext = i-1;
          }
      }
  if (m_verbosity>1) xout << "Next iteration slot "<<m_iNext<<std::endl;

  //     bool
  //         PrintDIISstate = true;
  //     if ( PrintDIISstate ) {
  //         std::ostream &xout = xout;
  //         xout << "  iUsedVecs: "; for ( uint i = 0; i < nDim; ++ i ) xout << " " << iUsedVecs[i]; xout << std::endl;
  //         PrintMatrixGen( xout, m_ErrorMatrix.data(), nMaxDim(), 1, nMaxDim(), m_ErrorMatrix.nStride(), "DIIS-B (resident)" );
  //         PrintMatrixGen( xout, B.data(), nDim+1, 1, nDim+1, B.nStride(), "DIIS-B/I" );
  //         PrintMatrixGen( xout, Rhs.data(), 1, 1, nDim+1, 1, "DIIS-Rhs" );
  //         PrintMatrixGen( xout, Coeffs.data(), 1, 1, nDim+1, 1, "DIIS-C" );
  //         xout << std::endl;
  //     }
  //  xout << "before extrapolation "<<std::endl;
  //  for (auto vs=m_solutions.begin(); vs!=m_solutions.end(); vs++)
  //      xout << "solution vector"<<*vs<<std::endl;
  //      xout << "front solution vector"<<m_solutions.front()<<std::endl;
  //      xout << "0 solution vector"<<m_solutions[0]<<std::endl;

  // now actually perform the extrapolation on the residuals
  // and amplitudes.
  m_LastAmplitudeCoeff = Coeffs[iThis];
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
  for (size_t l=0; l<other.size(); l++) other[l]->zero();
  for (size_t k=0; k<Coeffs.rows(); k++) {
      //      xout << "extrapolation k="<<k<<",iUsedVecs[k]="<<iUsedVecs[k]<<std::endl;
      //      xout << "add solution "<<iUsedVecs[k]<<m_solutions[iUsedVecs[k]].front()<<" with factor "<<Coeffs[k]<<std::endl;
      residual.front()->axpy(Coeffs[k],m_residuals[iUsedVecs[k]].front());
      solution.front()->axpy(Coeffs[k],m_solutions[iUsedVecs[k]].front());
      for (size_t l=0; l<other.size(); l++)
        other[l]->axpy(Coeffs[k],m_others[iUsedVecs[k]][l]);
    }
  if (m_verbosity>2) {
      xout << "DIIS.extrapolate() final extrapolated solution: "<<solution.front()<<std::endl;
      xout << "DIIS.extrapolate() final extrapolated residual: "<<residual.front()<<std::endl;
    }
}

void DIIS::FindUsefulVectors(uint *iUsedVecs, uint &nDim, double &fBaseScale, uint iThis)
{
  // remove lines from the system which correspond to vectors which are too bad
  // to be really useful for extrapolation purposes, but which might break
  // numerical stability if left in.
  double const
      fThrBadResidual = 1e12;
  double
      fBestResidualDot = m_ErrorMatrix(iThis,iThis),
      fWorstResidualDot = fBestResidualDot;
  assert(m_iVectorAge[iThis] < VecNotPresent);
  for ( uint i = 0; i < m_iVectorAge.size(); ++ i ) {
      if ( m_iVectorAge[i] >= VecNotPresent ) continue;
      fBestResidualDot = std::min( m_ErrorMatrix(i,i), fBestResidualDot );
    }
  nDim = 0;
  for ( uint i = 0; i < m_iVectorAge.size(); ++ i ){
      if ( i != iThis && m_iVectorAge[i] >= VecNotPresent )
        continue;
      if ( i != iThis && m_ErrorMatrix(i,i) > fBestResidualDot * fThrBadResidual) {
          m_iVectorAge[i] = VecNotPresent; // ignore this slot next time.
          continue;
        }
      fWorstResidualDot = std::max( m_ErrorMatrix(i,i), fWorstResidualDot );
      iUsedVecs[nDim] = i;
      ++ nDim;
    }

  fBaseScale = std::sqrt(fWorstResidualDot * fBestResidualDot);
  if ( fBaseScale <= 0. )
    fBaseScale = 1.;
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
      converged = d.iterate(g,x);
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
  d.solve(g,x);
  xout  << "Distance from solution = "<<std::sqrt(((*x.front())[0]-1)*((*x.front())[0]-1)+((*x.front())[1]-1)*((*x.front())[1]-1)) <<std::endl;

}

