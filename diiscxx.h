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


class diisStorage;
class diis
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
  diis(std::vector<size_t> lengths, size_t maxDim=6, double threshold=1e-3, DiisMode_type DiisMode=DIIS, size_t buffer_size=1024);
  ~diis();
  void setOptions(size_t maxDim=6, double threshold=1e-3, enum DiisMode_type DiisMode=DIIS);
  /*!
   * \brief discards previous iteration vectors, but does not clear records
   */
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
  std::vector<diisStorage*> store_; // files for each vector
  double threshold_;
  size_t maxDim_;
  unsigned int nDim_;
  size_t buffer_size_;
  //> 0xffff: no vector in this slot. Otherwise: number of iterations
  // the vector in this slot has already been inside the DIIS system.
  std::vector<uint> m_iVectorAge;
  uint m_iNext; //< next vector to be overwritten. nDim+1 if nDim < MaxDim_.
  // find vectors which are not considered too bad for extrapolation purposes.
  void FindUsefulVectors(uint *iUsedVecs, uint &nDimUsed, double &fBaseScale, uint iThis);
  void InterpolateFrom(diisStorage* store_, double* result,  Eigen::VectorXd Coeffs, size_t length);

  Eigen::MatrixXd m_ErrorMatrix;
  std::vector<double> m_Weights;

  // the following variables are kept for informative/displaying purposes
  double
      // dot(R,R) of last residual vector fed into this state.
      m_LastResidualNormSq,
      // coefficient the actual new vector got in the last DIIS step
      m_LastAmplitudeCoeff;

  void dump(std::vector<double *> vectors, unsigned int index);
  void load(std::vector<double*> vectors, unsigned int index);
};

/*!
 * \brief The diisStorage class provides auxiliary storage (usually on an external file) for the diis class
 */
class diisStorage
{
public:
  /*!
   * \brief diisStorage
   * \param lengthHint A provided estimate of the total amount of storage that is likely to be needed.
   */
  diisStorage(size_t lengthHint=0);
  /*!
   * \brief Query the total storage used.
   */
  size_t size();
  virtual void write(const double* buffer, size_t length, size_t address);
  virtual void read(double* buffer, size_t length, size_t address);
  ~diisStorage();
private:
  std::fstream dumpFile_;
protected:
  size_t size_; // internally store number of bytes used
};

#endif // DIISCXX_H
