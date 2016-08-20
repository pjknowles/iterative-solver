#ifndef DIISCXX_H
#define DIISCXX_H
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
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
   * \param AmpLength
   * \param ResLength
   * \param OtherLength
   * \param maxDim Maximum DIIS dimension allowed
   * \param threshold Residual threshold for inclusion of a vector in the DIIS state.
   * \param DiisMode
   */
  diis(size_t AmpLength, size_t ResLength, size_t OtherLength, size_t maxDim=6, double threshold=1e-3, DiisMode_type DiisMode=DIIS);
//  diis(size_t AmpLength, size_t maxDim=6, double threshold=1e-3, DiisMode_type DiisMode=DIIS);
  ~diis();
  void setOptions(size_t maxDim=6, double threshold=1e-3, enum DiisMode_type DiisMode=DIIS);
  /*!
   * \brief discards previous iteration vectors, but does not clear records
   */
  void Reset();
  /*!
   * \brief Introduce a new iteration vector, and perform extrapolation
   * \param Amp
   * \param Res
   * \param Other
   * \param weight
   */
  void extrapolate (double *Amp, double* Res=nullptr, double* Other=nullptr, double weight=1.0);

  double fLastResidual() const { return m_LastResidualNormSq; }
  double fLastCoeff() const { return m_LastAmplitudeCoeff; }
  unsigned int nLastDim() const {return std::count(m_iVectorAge.begin(),m_iVectorAge.end(),VecNotPresent);}
  unsigned int nNextVec() const { return m_iNext; }
  unsigned int nMaxDim() const { return maxDim_; }
private:
  typedef unsigned int uint;
  diis();
  size_t AmpLength_;
  size_t ResLength_;
  size_t OtherLength_;
  size_t TotalLength_;
  enum DiisMode_type DiisMode_;
  diisStorage* store_;
  double threshold_;
  size_t maxDim_;
  unsigned int nDim_;
  //> 0xffff: no vector in this slot. Otherwise: number of iterations
  // the vector in this slot has already been inside the DIIS system.
  std::vector<uint> m_iVectorAge;
  uint m_iNext; //< next vector to be overwritten. nDim+1 if nDim < MaxDim_.
  // find vectors which are not considered too bad for extrapolation purposes.
  void FindUsefulVectors(uint *iUsedVecs, uint &nDimUsed, double &fBaseScale, uint iThis);
  void InterpolateFrom(double* result, double fOwnCoeff, Eigen::VectorXd Coeffs, size_t length, size_t offset);

  Eigen::MatrixXd m_ErrorMatrix;
  std::vector<double> m_Weights;

  // the following variables are kept for informative/displaying purposes
  double
      // dot(R,R) of last residual vector fed into this state.
      m_LastResidualNormSq,
      // coefficient the actual new vector got in the last DIIS step
      m_LastAmplitudeCoeff;

  void dump(double *Amp, double* Res, double* Other, unsigned int index);
  void load(double *Amp, double* Res, double* Other, unsigned int index);
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
