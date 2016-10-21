#ifndef RSPT_H
#define RSPT_H
#include "IterativeSolver.h"

namespace IterativeSolver{

  /** @example RSPTexample.cpp */
  /*!
 * \brief A class that finds the lowest eigensolution of a matrix as a perturbation series
 *
 * Example of simplest use: @include RSPTexample.cpp
 *
 */
class RSPT : public IterativeSolverBase
{
  public:
    /*!
     * \brief RSPT
   * \param residualFunction A function that evaluates the residual vectors. Used by method solve(); does not have to be provided if iterations are constructed explicitly in the calling program.
   * \param preconditionerFunction A function that applies a preconditioner to a residual to give an update. Used by methods iterate() and solve().
   * Also used to construct the action of H0 if the shift parameter is zero.
     */
    RSPT(const ParameterSetTransformation residualFunction, const ParameterSetTransformation preconditionerFunction=&IterativeSolver::steepestDescent);
    static void test (size_t n, double alpha);
  protected:
    virtual void extrapolate(ParameterVectorSet & residual, ParameterVectorSet & solution, ParameterVectorSet & other, const optionMap options=optionMap());
  public:
    int m_order; ///< Up to what order of perturbation theory should the energy be obtained.
    std::vector<double> incremental_energies(size_t state=0) {return m_incremental_energies;} ///< The incremental energies order by order.
    double energy(size_t order, size_t state=0); ///< The total energy to a given order.
  private:
    RSPT();
    ParameterVectorSet m_lastH0mE0psi;
    std::vector<double> m_incremental_energies;
    double m_E0;
};

}


#endif // RSPT_H
