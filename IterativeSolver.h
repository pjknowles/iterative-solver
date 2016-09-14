#ifndef ITERATIVESOLVER_H
#define ITERATIVESOLVER_H
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
   * \brief Set convergence threshold
   * \param thresh If predicted residual . solution is less than this, converged, irrespective of cthresh and gthresh.
   */
  void setThresholds(double thresh) { thresh_=thresh;}
protected:
  double calculateError(const ParameterVectorSet & residual, const ParameterVectorSet & solution);
  std::vector<double> calculateErrors(const ParameterVectorSet & residual, const ParameterVectorSet & solution);
  int verbosity_;
  double thresh_;
  std::vector<ParameterVectorSet> residuals_;
  std::vector<ParameterVectorSet> solutions_;
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

#endif // ITERATIVESOLVER_H
