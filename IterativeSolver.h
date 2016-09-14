#ifndef DIISCXX_H
#define DIISCXX_H
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <numeric>
#include "ParameterVector.h"
#include <Eigen/Dense>

#ifndef nullptr
#define nullptr NULL
#endif
#define VecNotPresent 0xffff

namespace IterativeSolver {

class Storage;

/*!
 * \brief Place-holding template for residual calculation. It just returns the input as output.
 * \param inputs The parameters.
 * \param outputs On output, contains the corresponding residual vectors.
 * \param shift
 */
void noOp(const ParameterVectorSet & inputs, ParameterVectorSet & outputs, ParameterScalar shift=0) { outputs=inputs; }
/*!
 * \brief Place-holding template for update calculation. It just returns the input as output.
 * \param inputs The parameters.
 * \param outputs On output, contains the corresponding residual vectors.
 * \param shift
 */
void steepestDescent(const ParameterVectorSet & inputs, ParameterVectorSet & outputs, ParameterScalar shift=0) {
    for (size_t k=0; k<inputs.size(); k++)
        outputs[k].axpy(-1,inputs[k]);
}
/*!
 * \brief A base class for iterative solvers such as DIIS, KAIN, Davidson. The class provides support for preconditioned update, via a provided function.
 */
class IterativeSolverBase
{
protected:
  typedef void (*ParameterSetTransformation)(const ParameterVectorSet & inputs, ParameterVectorSet & outputs, ParameterScalar shift);
  /*!
   * \brief IterativeSolverBase
   * \param updateFunction A function that applies a preconditioner to a residual to give an update. Used by methods iterate() and solve().
   * \param residualFunction A function that evaluates the residual vectors. Used by method solve(); does not have to be provided if iterations are constructed explicitly in the calling program.
   */
  IterativeSolverBase(ParameterSetTransformation updateFunction=&IterativeSolver::steepestDescent, ParameterSetTransformation residualFunction=&IterativeSolver::noOp);
  virtual ~IterativeSolverBase();
protected:
  /*!
   * \brief The function that will take the current solution and residual, and produce the predicted solution.
   */
  ParameterSetTransformation updateFunction_;
  /*!
   * \brief The function that will take a current solution and calculate the residual.
   */
  ParameterSetTransformation residualFunction_;
public:
  /*!
   * \brief Take, typically, a current solution and residual, and return new solution.
   * In the context of Lanczos-like methods, the input will be a current expansion vector and the result of
   * acting on it with the matrix, and the output will be a new expansion vector.
   * \param residual
   * \param solution
   */
  virtual bool iterate(const ParameterVectorSet & residual, ParameterVectorSet & solution);
  /*!
   * \brief Solve iteratively by repeated calls to residualFunction() and iterate().
   * \param residual Ignored on input; on exit, contains the final residual; used as working space.
   * \param solution On input, contains an initial guess; on exit, contains the final solution.
   * \return Whether or not convergence has been reached.
   */
  void solve(ParameterVectorSet & residual, ParameterVectorSet & solution);
  /*!
   * \brief Control level of output
   * \param verbosity The higher the number, the more output.
   */
  void setVerbosity(int verbosity) { verbosity_=verbosity;}
  /*!
   * \brief Set convergence thresholds
   * \param ethresh If predicted residual . solution is less than this, converged, irrespective of cthresh and gthresh.
   * \param cthresh If norm of update is less than this, converged, irrespective of ethresh and gthresh.
   * \param gthresh If norm of residual is less than this, converged, irrespective ofecthresh and cthresh.
   */
  void setThresholds(double ethresh, double cthresh=0, double gthresh=0) { ethresh_=ethresh; cthresh_=cthresh; gthresh_=gthresh;}
protected:
  int verbosity_;
  double ethresh_;
  double cthresh_;
  double gthresh_;
};

/*!
 * \brief A class that encapsulates accelerated convergence of non-linear equations through the DIIS or related methods.
 *
 * Example of simplest use, with DIIS extrapolation based on the residual as error vector:
 * \code
 * std::vector<double> x(2);
 * std::vector<double> g(2);
 * x[0]=x[1]=0.9; // initial guess
 * std::vector<double> diag(2); diag[0]=700; diag[1]=200; // preconditioner
 * std::vector<size_t> lengths; lengths.push_back(g.size()); lengths.push_back(x.size());
 * Diis d(lengths);
 * d.addPreconditioner(&diag[0],0,true);
 * for (int iteration=1; iteration < 1000 && d.fLastResidual() > 1e-25; iteration++) {
 *     g[0]=2*x[0]-2+400*x[0]*(x[0]*x[0]-x[1]); g[1]=200*(x[1]-x[0]*x[0]); // Rosenbrock function gradient
 *     d.iterate(&g[0],&x[0]);
 *     std::cout << "iteration "<<iteration<<", Residual norm = "<<std::sqrt(d.fLastResidual())
 *                 << ", Distance from solution = "<<std::sqrt((x[0]-1)*(x[0]-1)+(x[1]-1)*(x[1]-1))<<std::endl;
 *   }
 * \endcode
 */
class Diis : public IterativeSolverBase
{
public:
  enum DiisMode_type {disabled ///< No extrapolation is performed
                      , DIIS ///< Direct Inversion in the Iterative Subspace
                      , KAIN ///< Krylov Accelerated Inexact Newton
                     };
/*!
   * \brief diis
   * \param lengths A list of lengths of vectors that will be extrapolated.
   * The first vector
   * will be the one that is analysed to construct the extrapolation, and the remainder only have
   * the extrapolation applied to them.
   * \param maxDim Maximum DIIS dimension allowed
   * \param threshold Residual threshold for inclusion of a vector in the DIIS state.
   * \param DiisMode Whether to perform DIIS, KAIN, or nothing.
   * \param globalSum A function that in a parallel environment sums globally, then broadcasts, the vector supplied to it.
   */
  Diis(ParameterSetTransformation updateFunction=&IterativeSolver::steepestDescent, ParameterSetTransformation residualFunction=&IterativeSolver::noOp);
  ~Diis();
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
   * \param weight Weight to be given to this vector.
   */
  void extrapolate (ParameterVectorSet vectors, double weight=1.0);

  /*!
   * \brief Perform DIIS extrapolation based on a residual vector, and then update the extrapolated solution with the extrapolated residual
   * \param residual The residual vector. On exit, it contains the extrapolated residual.
   * \param solution The current solution vector. On exit, it contains the extrapolated and updated solution.
   * \param weight Weight to be given to this vector.
   * \param other Any other vectors that should be extrapolated.
   *
   * The function performs extrapolate() followed by update().
   * The lengths of the vectors must have been declared in the object constructor,
   * with the two (equal) lengths of residual and solution coming before, in order, any contents of other.
   */
  void iterate (double* residual, double* solution, double weight=1.0, std::vector<double*> other=std::vector<double*>(0));

  /*!
   * \brief Return the square L2 norm of the extrapolated residual from the last call to extrapolate() or iterate().
   * \return
   */
  double fLastResidual() const { return m_LastResidualNormSq; }
  /*!
   * \brief Return the coefficient of the last residual vector in the extrapolated residual from the last call to extrapolate() or iterate().
   * \return
   */
  double fLastCoeff() const { return m_LastAmplitudeCoeff; }
  unsigned int nLastDim() const {return std::count(m_iVectorAge.begin(),m_iVectorAge.end(),VecNotPresent);}
  unsigned int nNextVec() const { return m_iNext; }
  unsigned int nMaxDim() const { return maxDim_; }
private:
  typedef unsigned int uint;
  Diis();
  std::vector<size_t> lengths_;
  size_t totalLength() { return std::accumulate(lengths_.begin(),lengths_.end(),0); }
  enum DiisMode_type DiisMode_;
  std::vector<Storage*> store_; // files for each vector
  double threshold_;
  size_t maxDim_;
  unsigned int nDim_;
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
