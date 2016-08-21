#include "IterativeSolver.h"
#include <algorithm>

using namespace IterativeSolver;

IterativeSolverBase::IterativeSolverBase() : buffer_size_(1024) {}

diis::diis(std::vector<size_t> lengths, size_t maxDim, double threshold, DiisMode_type DiisMode, size_t buffer_size)
  : lengths_(lengths), maxDim_ (maxDim), threshold_(threshold), DiisMode_(DiisMode), verbosity_(-1)
{
  buffer_size_=buffer_size;
  setOptions(maxDim,threshold,DiisMode);
  for (std::vector<size_t>::const_iterator l=lengths.begin(); l!=lengths.end(); l++)
    store_.push_back(new Storage(maxDim*(*l)*sizeof(double)));
  m_LastResidualNormSq=1e99; // so it can be tested even before extrapolation is done
  Reset();
}

diis::~diis()
{
  for (std::vector<Storage*>::const_iterator s=store_.begin(); s!=store_.end(); s++)
    delete *s;
}

void diis::setOptions(size_t maxDim, double threshold, DiisMode_type DiisMode)
{
   maxDim_ = maxDim;
   threshold_ = threshold;
   DiisMode_ = DiisMode;
   Reset();
//   std::cout << "maxDim_ set to "<<maxDim_<<" in setOptions"<<std::endl;
}

void diis::Reset()
{

}

double Dot(double* a, double* b, size_t n)
{
  double result=0;
  for (size_t k=0; k<n; k++) result+=a[k]*b[k];
  return result;
}

void diis::LinearSolveSymSvd(Eigen::VectorXd& Out, const Eigen::MatrixXd& Mat, const Eigen::VectorXd& In, unsigned int nDim, double Thr)
{
    Eigen::VectorXd Ews(nDim);
    Eigen::VectorXd Xv(nDim); // input vectors in EV basis.
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    es.compute(-Mat);
    Ews = es.eigenvalues();

    if (verbosity_ > 1) {
        std::cout << "diis::LinearSolveSymSvd"<<std::endl;
        std::cout << "Mat="<<Mat<<std::endl;
        std::cout << "Ews="<<Ews<<std::endl;
        std::cout << "In="<<In<<std::endl;
        std::cout << "es.eigenvectors()="<<es.eigenvectors()<<std::endl;
      }

    Xv = In.transpose()*es.eigenvectors();
//    std::cout << "Xv="<<Xv<<std::endl;
    for (size_t iEw = 0; iEw != nDim; ++ iEw)
        if (std::abs(Ews(iEw)) >= Thr)
            Xv(iEw) /= -Ews(iEw);
            // ^- note that this only screens by absolute value.
            // no positive semi-definiteness is assumed!
        else
            Xv(iEw) = 0.;
    if (verbosity_ > 1) std::cout << "Xv="<<Xv<<std::endl;
    Out = es.eigenvectors() * Xv;
    if (verbosity_ > 1) std::cout << "Out="<<Out<<std::endl;
    Out=Out.block(0,0,Out.size()-1,1);
//    std::cout << "Out="<<Out<<std::endl;
}



void diis::extrapolate (std::vector<double*> vectors, double weight)
{
  if (maxDim_ <= 1 || DiisMode_ == disabled) return;

  double fThisResidualDot = Dot(vectors[0],vectors[0],lengths_[0]);
  m_LastResidualNormSq = fThisResidualDot;

  if ( m_iNext == 0 && fThisResidualDot > threshold_ )
      // current vector is to be considered too wrong to be useful for DIIS
      // purposes. Don't store it.
      return;

  uint iThis = m_iNext;
  if (verbosity_ > 1) std::cout<< "iThis=m_iNext "<<m_iNext<<std::endl;
  assert(iThis < maxDim_);
  if (iThis >= m_iVectorAge.size()) m_iVectorAge.resize(iThis+1);
  if (iThis >= m_ErrorMatrix.cols()) m_ErrorMatrix.conservativeResize(iThis+1,iThis+1);
  m_ErrorMatrix(iThis,iThis)=fThisResidualDot;
  for ( uint i = 0; i < m_iVectorAge.size(); ++ i )
      m_iVectorAge[i] += 1;
  m_iVectorAge[iThis] = 0;
  if (verbosity_>1) {
      std::cout << "iVectorAge:";
      for (std::vector<uint>::const_iterator a=m_iVectorAge.begin(); a!=m_iVectorAge.end(); a++)
        std::cout << " "<<*a;
      std::cout << std::endl;
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
  if (verbosity_>1) {
      std::cout << "iUsedVecs:";
      for (std::vector<uint>::const_iterator a=iUsedVecs.begin(); a!=iUsedVecs.end(); a++)
        std::cout << " "<<*a;
      std::cout << std::endl;
    }
  // transform iThis into a relative index.
  for ( uint i = 0; i < nDim; ++ i )
      if ( iThis == iUsedVecs[i] ) {
          iThis = i;
          break;
      }

  // write current residual and other vectors to their designated place
  if (verbosity_>0) std::cout << "write current vectors to record "<<iUsedVecs[iThis]<<std::endl;
  for (unsigned int k=0; k<lengths_.size(); k++)
   store_[k]->write(vectors[k],lengths_[k]*sizeof(double),iUsedVecs[iThis]*lengths_[k]*sizeof(double));
  if (m_Weights.size()<=iUsedVecs[iThis]) m_Weights.resize(iUsedVecs[iThis]+1);
  m_Weights[iUsedVecs[iThis]] = weight;

  // go through previous residual vectors and form the dot products with them
//  std::vector<double> ResDot(nDim);
//  {
//  std::vector<double> buffer(std::min(buffer_size_,lengths_[0]));
//  for (uint i=0; i<nDim; i++) {
//      ResDot[i] = 0;
//      for (uint block=0; block<lengths_[0]; block+=buffer.size()) {
//          store_[0]->read(&buffer[0],std::min(buffer.size(),lengths_[0]-block)*sizeof(double),(i*lengths_[0]+block)*sizeof(double));
//          ResDot[i] += Dot(&vectors[0][block],&buffer[0],std::min(buffer.size(),lengths_[0]-block));
//        }
//    }
//  }
  std::vector<double> ResDot = overlapsWithStore(store_[0], vectors[0], lengths_[0], nDim);
  if (iThis >= nDim) ResDot.resize(iThis+1);
  ResDot[iThis] = fThisResidualDot;


  // update resident error matrix with new residual-dots
  for ( uint i = 0; i < nDim; ++ i ) {
      m_ErrorMatrix(iUsedVecs[i], iUsedVecs[iThis]) = ResDot[i];
      m_ErrorMatrix(iUsedVecs[iThis], iUsedVecs[i]) = ResDot[i];
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
  if (verbosity_>1) std::cout << "Combination of iteration vectors: "<<Coeffs.transpose()<<std::endl;

  // Find a storage place for the vector in the next round. Either
  // an empty slot or the oldest vector.
  uint iOldestAge = m_iVectorAge[0];
  m_iNext = 0;
  if (m_iVectorAge.size() < maxDim_) {
    m_iNext = m_iVectorAge.size(); m_iVectorAge.push_back(0);
    } else
    for ( uint i = m_iVectorAge.size(); i != 0; -- i ){
      if ( iOldestAge <= m_iVectorAge[i-1] ) {
          iOldestAge = m_iVectorAge[i-1];
          m_iNext = i-1;
      }
  }

//     bool
//         PrintDiisState = true;
//     if ( PrintDiisState ) {
//         std::ostream &xout = std::cout;
//         xout << "  iUsedVecs: "; for ( uint i = 0; i < nDim; ++ i ) xout << " " << iUsedVecs[i]; xout << std::endl;
//         PrintMatrixGen( xout, m_ErrorMatrix.data(), nMaxDim(), 1, nMaxDim(), m_ErrorMatrix.nStride(), "DIIS-B (resident)" );
//         PrintMatrixGen( xout, B.data(), nDim+1, 1, nDim+1, B.nStride(), "DIIS-B/I" );
//         PrintMatrixGen( xout, Rhs.data(), 1, 1, nDim+1, 1, "DIIS-Rhs" );
//         PrintMatrixGen( xout, Coeffs.data(), 1, 1, nDim+1, 1, "DIIS-C" );
//         xout << std::endl;
//     }

  // now actually perform the extrapolation on the residuals
  // and amplitudes.
  m_LastAmplitudeCoeff = Coeffs[iThis];
//  TR.InterpolateFrom( Coeffs[iThis], Coeffs.data(), ResRecs, AmpRecs,
//      nDim, iThis, *m_Storage.pDevice, m_Memory );
  for (uint k=0; k<vectors.size(); k++) {
   InterpolateFrom(store_[k],vectors[k],Coeffs,lengths_[k],iUsedVecs);
}
}

std::vector<double> IterativeSolverBase::overlapsWithStore(Storage* store, double* vector, size_t length , size_t nVector)
{
  std::vector<double> result(nVector);
  std::vector<double> buffer(std::min(buffer_size_,length));
  for (unsigned int i=0; i<nVector; i++) {
      result[i] = 0;
      for (unsigned int block=0; block<length; block+=buffer.size()) {
          store->read(&buffer[0],std::min(buffer.size(),length-block)*sizeof(double),(i*length+block)*sizeof(double));
          result[i] += Dot(&vector[block],&buffer[0],std::min(buffer.size(),length-block));
        }
    }
  return result;
}
void IterativeSolverBase::InterpolateFrom(Storage* store_, double* result,  Eigen::VectorXd Coeffs, size_t length, std::vector<unsigned int>& iUsedVecs)
{
  std::vector<double> buffer(std::min(buffer_size_,length));
  for (size_t i=0; i<length; i++) result[i]=0;
  for (unsigned int i=0; i<Coeffs.rows(); i++) {
      for (unsigned int block=0; block<length; block+=buffer.size()) {
          store_->read(&buffer[0],std::min(buffer.size(),length-block)*sizeof(double),(iUsedVecs[i]*length+block)*sizeof(double));
          for (size_t j=0; j<std::min(buffer.size(),length-block); j++) {
              result[j+block] += buffer[j] * Coeffs[i];
            }
        }
    }
}

void diis::FindUsefulVectors(uint *iUsedVecs, uint &nDim, double &fBaseScale, uint iThis)
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


// default I/O implementation
#include <iostream>
#include <cstdio>
#include <cstring>
#include <cstdlib>
Storage::Storage(size_t lengthHint, int option)
{
  char *tmpname=strdup("tmpfileXXXXXX");
  mkstemp(tmpname);
  dumpFile_.open (tmpname, std::ios::out | std::ios::in | std::ios::binary);
  free(tmpname);
  size_=0;
}

Storage::~Storage()
{
  dumpFile_.close();
}

void Storage::write(const double *buffer, size_t length, size_t address)
{
  dumpFile_.seekg(address);//,std::ios::beg);
  dumpFile_.write((char*) buffer,length);
//  std::cout << "Storage::write at address "<< address << ", length="<<length<<","<< buffer[0]<<std::endl;
  if (length+address > size_) size_ = length+address;
}

void Storage::read(double *buffer, size_t length, size_t address)
{
  if (address+length > size_) throw std::range_error("Storage: attempt to load from beyond end of storage");
  dumpFile_.seekg(address);//,std::ios::beg);
  dumpFile_.read((char*) buffer,length);
//  std::cout << "Storage::read at address "<< address << ", length="<<length<<","<< buffer[0]<<std::endl;
}

size_t Storage::size()
{
  return size_;
}
