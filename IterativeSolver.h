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

/*!
 * \brief Place-holding template for residual calculation. It just returns the input as output.
 * \param inputs The parameters.
 * \param outputs On output, contains the corresponding residual vectors.
 * \param shift
 */
void noOp(const ParameterVectorSet & inputs, ParameterVectorSet & outputs, std::vector<ParameterScalar> shift=std::vector<ParameterScalar>()) { outputs=inputs; }
/*!
 * \brief Place-holding template for update calculation. It just returns the input as output.
 * \param inputs The parameters.
 * \param outputs On output, contains the corresponding residual vectors.
 * \param shift
 */
void steepestDescent(const ParameterVectorSet & inputs, ParameterVectorSet & outputs, std::vector<ParameterScalar> shift=std::vector<ParameterScalar>()) {
    for (size_t k=0; k<inputs.size(); k++)
        outputs[k].axpy(-1,inputs[k]);
}
/*!
 * \brief A base class for iterative solvers such as DIIS, KAIN, Davidson. The class provides support for preconditioned update, via a provided function.
 */
class IterativeSolverBase
{
protected:
  typedef void (*ParameterSetTransformation)(const ParameterVectorSet & inputs, ParameterVectorSet & outputs, std::vector<ParameterScalar> shift);
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
  ParameterSetTransformation m_updateFunction;
  /*!
   * \brief The function that will take a current solution and calculate the residual.
   */
  ParameterSetTransformation m_residualFunction;
public:
  /*!
   * \brief Take, typically, a current solution and residual, and return new solution.
   * In the context of Lanczos-like methods, the input will be a current expansion vector and the result of
   * acting on it with the matrix, and the output will be a new expansion vector.
   * iterate() saves the vectors, calls extrapolate(), calls m_updateFunction(), calls calculateErrors(), and then assesses the error.
   * Derivative classes may often be able to be implemented by changing only extrapolate(), not iterate() or solve().
   * \param residual On input, the residual for solution on entry. On exit, the extrapolated residual.
   * \param solution On input, the current solution or expansion vector. On exit, the next solution or expansion vector.
   * \param other Optional additional vectors that should be extrapolated.
   * \param options A string of options to be interpreted by extrapolate().
   */
  virtual bool iterate(ParameterVectorSet & residual, ParameterVectorSet & solution, ParameterVectorSet & other, std::string options="");
  virtual bool iterate(ParameterVectorSet & residual, ParameterVectorSet & solution, std::string options="") { ParameterVectorSet other; return iterate(residual,solution,other,options); }
  /*!
   * \brief Solve iteratively by repeated calls to residualFunction() and iterate().
   * \param residual Ignored on input; on exit, contains the final residual; used as working space.
   * \param solution On input, contains an initial guess; on exit, contains the final solution.
   * \return Whether or not convergence has been reached.
   */
  virtual void solve(ParameterVectorSet & residual, ParameterVectorSet & solution);
  /*!
   * \brief Control level of output
   * \param verbosity The higher the number, the more output.
   */
  void setVerbosity(int verbosity) { m_verbosity=verbosity;}
  /*!
   * \brief Set convergence threshold
   * \param thresh If predicted residual . solution is less than this, converged, irrespective of cthresh and gthresh.
   */
  void setThresholds(double thresh) { m_thresh=thresh;}
protected:
  virtual void extrapolate(ParameterVectorSet & residual, ParameterVectorSet & solution, ParameterVectorSet & other, std::string options="");
  virtual void extrapolate(ParameterVectorSet & residual, ParameterVectorSet & solution, std::string options="") { ParameterVectorSet other; extrapolate(residual,solution,other,options); }
  double calculateError(const ParameterVectorSet & residual, const ParameterVectorSet & solution);
  std::vector<double> calculateErrors(const ParameterVectorSet & residual, const ParameterVectorSet & solution);
  int m_verbosity;
  double m_thresh;
  std::vector<ParameterVectorSet> m_residuals;
  std::vector<ParameterVectorSet> m_solutions;
  std::vector<ParameterVectorSet> m_others;
};
}


#endif // ITERATIVESOLVER_H
