#ifndef DIIS_H
#define DIIS_H

#include "IterativeSolver.h"

namespace IterativeSolver {

#define VecNotPresent 0xffff

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
   */
  Diis(ParameterSetTransformation updateFunction=&IterativeSolver::steepestDescent, ParameterSetTransformation residualFunction=&IterativeSolver::noOp);
  ~Diis();
  /*!
   * \brief Set options for DIIS.
   * \param maxDim Maximum DIIS dimension allowed
   * \param threshold Residual threshold for inclusion of a vector in the DIIS state.
   * \param DiisMode Whether to perform DIIS, KAIN, or nothing.
   */
  void setOptions(size_t maxDim=6, double threshold=1e-3, enum DiisMode_type DiisMode=DIIS);
  /*!
   * \brief discards previous iteration vectors, but does not clear records
   */
  void Reset();
  /*!
   * \brief Introduce a new iteration vector, and perform extrapolation
   * \param residual
   * The vector that
   * will be the one that is analysed to construct the extrapolation.
   * \param solution
   * The current solution that gave rise to residual, and which will be extrapolated to a new predicted solution.
   * \param other (optional)
   * Corresponding other vectors whose sequence will be extrapolated.
   * \param options can contain "weight=xxx" where xxx is the weight to be given to this vector.
   */
  void extrapolate(ParameterVectorSet & residual, ParameterVectorSet & solution, ParameterVectorSet & other, std::string options="");

  /*!
   * \brief Perform DIIS extrapolation based on a residual vector, and then update the extrapolated solution with the extrapolated residual
   * The function performs extrapolate() followed by update().
   * \brief Introduce a new iteration vector, and perform extrapolation
   * \param residual
   * The vector that
   * will be the one that is analysed to construct the extrapolation.
   * \param solution
   * The current solution that gave rise to residual, and which will be extrapolated and updated to a new predicted solution.
   * \param other (optional)
   * Corresponding other vectors whose sequence will be extrapolated.
   * \param options can contain "weight=xxx" where xxx is the weight to be given to this vector.
   */
//  bool iterate(ParameterVectorSet & residual, ParameterVectorSet & solution, ParameterVectorSet & other, std::string options="");
//  bool iterate(ParameterVectorSet & residual, ParameterVectorSet & solution, std::string options="") { ParameterVectorSet other; return iterate(residual,solution,other,options); } // FIXME why necessary? Inherit?

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
  enum DiisMode_type DiisMode_;
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

}

#endif // DIIS_H