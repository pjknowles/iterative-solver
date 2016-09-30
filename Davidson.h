#ifndef DAVIDSON_H
#define DAVIDSON_H
#include "IterativeSolver.h"


namespace IterativeSolver{

  class Davidson : public IterativeSolverBase
  {
  public:
    Davidson();
    Davidson(ParameterSetTransformation updateFunction=&IterativeSolver::steepestDescent, ParameterSetTransformation residualFunction=&IterativeSolver::noOp);
    static void test(size_t dimension, size_t roots=1, int verbosity=0, int problem=0);
  protected:
    virtual void extrapolate(ParameterVectorSet & residual, ParameterVectorSet & solution, ParameterVectorSet & other, const optionMap options=optionMap());
  public:
    int m_roots; ///< How many roots to calculate (defaults to size of solution and residual vectors)
    std::vector<double> eigenvalues(); ///< The calculated eigenvalues
  };
}

#endif // DAVIDSON_H
