#ifndef DAVIDSON_H
#define DAVIDSON_H
#include "IterativeSolver.h"


namespace IterativeSolver{

  /** @example DavidsonExample.cpp */
  /*!
 * \brief A class that finds the lowest eigensolutions of a matrix using Davidson's method
 *
 * Example of simplest use: @include DavidsonExample.cpp
 *
 */
  class Davidson : public IterativeSolverBase
  {
  public:
    Davidson(const ParameterSetTransformation residualFunction, const ParameterSetTransformation preconditionerFunction=&IterativeSolver::steepestDescent);
    static void test(size_t dimension, size_t roots=1, int verbosity=0, int problem=0, bool orthogonalize=true);
  protected:
    virtual void extrapolate(ParameterVectorSet & residual, ParameterVectorSet & solution, ParameterVectorSet & other, const optionMap options=optionMap());
    virtual void report();
    Davidson();
  };
}

#endif // DAVIDSON_H