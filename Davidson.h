#ifndef DAVIDSON_H
#define DAVIDSON_H
#include "IterativeSolver.h"


namespace IterativeSolver{

class Davidson : public IterativeSolverBase
{
public:
    Davidson();
    Davidson(ParameterSetTransformation updateFunction=&IterativeSolver::steepestDescent, ParameterSetTransformation residualFunction=&IterativeSolver::noOp);
    public:
      /*!
       * \brief Take, a current expansion vector and residual, and return new expansion vector.
       * \param residual
       * \param solution
       * \return Whether or not convergence has been reached.
       */
      bool iterate(const ParameterVectorSet & residual, ParameterVectorSet & solution);
      /*!
       * \brief Solve iteratively by repeated calls to residualFunction() and iterate().
       * \param residual Ignored on input; on exit, contains the final residual; used as working space.
       * \param solution On input, contains an initial guess; on exit, contains the final solution.
       */
      void solve(ParameterVectorSet & residual, ParameterVectorSet & solution);
private:
    Eigen::MatrixXd m_SubspaceMatrix;
};
}

#endif // DAVIDSON_H
