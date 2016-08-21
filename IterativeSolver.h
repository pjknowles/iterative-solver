#ifndef DIISCXX_H
#define DIISCXX_H
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <Eigen/Dense>

#define nullptr NULL
#define VecNotPresent 0xffff

namespace IterativeSolver {

class Storage;
class IterativeSolverBase
{
protected:
  IterativeSolverBase();
  size_t buffer_size_;
  void StorageCombine(Storage* store_, double* result,  Eigen::VectorXd Coeffs, size_t length, std::vector<unsigned int> &iUsedVecs);
  std::vector<double> StorageDot(Storage* store, double* vector, size_t length , size_t nVector);
};
class diis : public IterativeSolverBase
{
public:
  enum DiisMode_type {disabled, DIIS, KAIN};
/*!
   * \brief diis
   * \param lengths A list of lengths of vectors that will be extrapolated.
   * The first vector
   * will be the one that is analysed to construct the extrapolation, and the remainder only have
   * the extrapolation applied to them.
   * \param maxDim Maximum DIIS dimension allowed
   * \param threshold Residual threshold for inclusion of a vector in the DIIS state.
   * \param DiisMode Whether to perform DIIS, KAIN, or nothing.
   * \param buffer_size The size of blocks used in accessing old vectors on backing store.
   */
  diis(std::vector<size_t> lengths, size_t maxDim=6, double threshold=1e6, DiisMode_type DiisMode=DIIS, size_t buffer_size=1024);
  ~diis();
  void setOptions(size_t maxDim=6, double threshold=1e-3, enum DiisMode_type DiisMode=DIIS);
  /*!
   * \brief discards previous iteration vectors, but does not clear records
   */
  void setVerbosity(int verbosity) { verbosity_=verbosity;}
  void Reset();
  /*!
   * \brief Introduce a new iteration vector, and perform extrapolation
   * \param vectors A list of vectors that will be extrapolated.
   * The first vector
   * will be the one that is analysed to construct the extrapolation, and the remainder only have
   * the extrapolation applied to them.
   * \param weight
   */
  void extrapolate (std::vector<double*> vectors, double weight=1.0);

  double fLastResidual() const { return m_LastResidualNormSq; }
  double fLastCoeff() const { return m_LastAmplitudeCoeff; }
  unsigned int nLastDim() const {return std::count(m_iVectorAge.begin(),m_iVectorAge.end(),VecNotPresent);}
  unsigned int nNextVec() const { return m_iNext; }
  unsigned int nMaxDim() const { return maxDim_; }
private:
  typedef unsigned int uint;
  diis();
  std::vector<size_t> lengths_;
  size_t totalLength() { std::accumulate(lengths_.begin(),lengths_.end(),0); }
  enum DiisMode_type DiisMode_;
  std::vector<Storage*> store_; // files for each vector
  double threshold_;
  size_t maxDim_;
  unsigned int nDim_;
  int verbosity_;
  //> 0xffff: no vector in this slot. Otherwise: number of iterations
  // the vector in this slot has already been inside the DIIS system.
  std::vector<uint> m_iVectorAge;
  uint m_iNext; //< next vector to be overwritten. nDim+1 if nDim < MaxDim_.
  // find vectors which are not considered too bad for extrapolation purposes.
  void FindUsefulVectors(uint *iUsedVecs, uint &nDimUsed, double &fBaseScale, uint iThis);
  void LinearSolveSymSvd(Eigen::VectorXd& Out, const Eigen::MatrixXd& Mat, const Eigen::VectorXd& In, unsigned int nDim, double Thr);

  Eigen::MatrixXd m_ErrorMatrix;
  std::vector<double> m_Weights;

  // the following variables are kept for informative/displaying purposes
  double
      // dot(R,R) of last residual vector fed into this state.
      m_LastResidualNormSq,
      // coefficient the actual new vector got in the last DIIS step
      m_LastAmplitudeCoeff;

};

/*!
 * \brief The Storage class provides auxiliary storage (usually on an external file) for
 * iterative solvers.
 */
class Storage
{
public:
  /*!
   * \brief Storage
   * \param lengthHint A provided estimate of the total amount of storage that is likely to be needed.
   * \param option An implementation-dependent parameter that controls operation
   */
  Storage(size_t lengthHint=0, int option=0);
  ~Storage();
  /*!
   * \brief Write data to the store.
   * \param buffer Provides the data to be written.
   * \param length Length of data, in bytes.
   * \param address Offset in store, in bytes.
   */
  virtual void write(const double* buffer, size_t length, size_t address);
  /*!
   * \brief Read data from the store.
   * \param buffer Receives the data to be read.
   * \param length Length of data, in bytes.
   * \param address Offset in store, in bytes.
   */
  virtual void read(double* buffer, size_t length, size_t address);
  /*!
   * \brief Query the total storage used.
   */
  virtual size_t size();
private:
  std::fstream dumpFile_;
  size_t size_; //< total storage (bytes) used
};

}

#endif // DIISCXX_H
