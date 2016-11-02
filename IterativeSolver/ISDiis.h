#ifndef DIIS_H
#define DIIS_H

#include "IterativeSolver.h"

namespace LinearAlgebra {

  /** @example DIISexample.cpp */
  /*!
 * \brief A class that encapsulates accelerated convergence of non-linear equations through the DIIS or related methods.
 *
 * Example of simplest use: @include DIISexample.cpp
 *
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
    DIIS(const ParameterSetTransformation residualFunction, const ParameterSetTransformation updateFunction=&LinearAlgebra::steepestDescent);
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
    void extrapolate(ParameterVectorSet & residual, ParameterVectorSet & solution, const optionMap options=optionMap()) { ParameterVectorSet other; extrapolate(residual,solution,other,options); }

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
    static void randomTest(size_t sample, size_t n=100, double alpha=0.1, double gamma=0.0, DIISmode_type mode=DIISmode);
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
    KAIN(const ParameterSetTransformation residualFunction=&LinearAlgebra::noOp, const ParameterSetTransformation updateFunction=&LinearAlgebra::steepestDescent)
      : DIIS(residualFunction,updateFunction) { setMode(KAINmode);}
//    void setMode( enum DIISmode_type mode=KAINmode)
//    { DIIS::setMode(mode); }
  private:
    KAIN();
  };

}

#endif // DIIS_H
