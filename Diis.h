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
 *     xout << "iteration "<<iteration<<", Residual norm = "<<std::sqrt(d.fLastResidual())
 *                 << ", Distance from solution = "<<std::sqrt((x[0]-1)*(x[0]-1)+(x[1]-1)*(x[1]-1))<<std::endl;
 *   }
 * \endcode
 */
  class DIIS : public IterativeSolverBase
  {
  public:
    enum DIISmode_type {disabled ///< No extrapolation is performed
                        , DIISmode ///< Direct Inversion in the Iterative Subspace
                        , KAINmode ///< Krylov Accelerated Inexact Newton
                       };
    /*!
   */
    DIIS(ParameterSetTransformation updateFunction=&IterativeSolver::steepestDescent, ParameterSetTransformation residualFunction=&IterativeSolver::noOp);
    ~DIIS();
    /*!
   * \brief Set options for DIIS.
   * \param mode Whether to perform DIIS, KAIN, or nothing.
   */
    virtual void setMode(enum DIISmode_type mode=DIISmode);
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
   * \param options can contain "weight=xxx" where xxx is the weight to be given to this vector. These options would normally be passed as the corresponding parameter in iterate().
   */
  protected:
    void extrapolate(ParameterVectorSet & residual, ParameterVectorSet & solution, ParameterVectorSet & other, const optionMap options=optionMap());

    /*!
   * \brief Return the square L2 norm of the extrapolated residual from the last call to extrapolate() or iterate().
   * \return
   */
  public:
    double fLastResidual() const { return m_LastResidualNormSq; }
    /*!
   * \brief Return the coefficient of the last residual vector in the extrapolated residual from the last call to extrapolate() or iterate().
   * \return
   */
    double fLastCoeff() const { return m_LastAmplitudeCoeff; }
    unsigned int nLastDim() const {return m_residuals.size();}
    unsigned int nMaxDim() const { return m_maxDim; }
    /*!
   * \brief Test the correct operation of the module. If an error is found, an exception is thrown.
   * \param verbosity How much to print.
   * - -1 Nothing at all is printed.
   * - 0 (default) Just a message that the test is taking place.
   * - 1, 2, 3,... more detail.
   * \param maxDim Maximum DIIS dimension allowed
   * \param svdThreshold Residual threshold for inclusion of a vector in the DIIS state.
   * \param mode Whether to perform DIIS, KAIN, or nothing.
   * \param difficulty Level of numerical challenge, ranging from 0 to 1.
   */
    static void test(int verbosity=0,
                     size_t maxDim=6,
                     double svdThreshold=1e-10,
                     enum DIISmode_type mode=DIISmode,
                     double difficulty=0.1);
    double m_svdThreshold; ///< Threshold for singular-value truncation in linear equation solver.
    size_t m_maxDim; ///< Maximum DIIS dimension allowed.
  private:
    typedef unsigned int uint;
    DIIS();
    enum DIISmode_type m_DIISmode;

    std::vector<double> m_Weights;

    // the following variables are kept for informative/displaying purposes
    double
    // dot(R,R) of last residual vector fed into this state.
    m_LastResidualNormSq,
    // coefficient the actual new vector got in the last DIIS step
    m_LastAmplitudeCoeff;

  };

  class KAIN : public DIIS
  {
  public:
    KAIN(ParameterSetTransformation updateFunction=&IterativeSolver::steepestDescent, ParameterSetTransformation residualFunction=&IterativeSolver::noOp)
      : DIIS(updateFunction,residualFunction) {}
    void setMode( enum DIISmode_type mode=KAINmode)
    { DIIS::setMode(mode); }
  private:
    KAIN();
  };

}

#endif // DIIS_H
