#ifndef DAVIDSON_H
#define DAVIDSON_H
#include "IterativeSolver.h"


namespace IterativeSolver{

class Davidson : public IterativeSolver
{
public:
    Davidson();
    public:
      /*!
       * \brief Take, typically, a current solution and residual, and return new solution.
       * In the context of Lanczos-like methods, the input will be a current expansion vector and the result of
       * acting on it with the matrix, and the output will be a new expansion vector.
       * \param residual
       * \param solution
       */
      bool iterate(const ParameterVectorSet & residual, ParameterVectorSet & solution);
      /*!
       * \brief Solve iteratively by repeated calls to residualFunction() and iterate().
       * \param residual Ignored on input; on exit, contains the final residual; used as working space.
       * \param solution On input, contains an initial guess; on exit, contains the final solution.
       * \return Whether or not convergence has been reached.
       */
      void solve(ParameterVectorSet & residual, ParameterVectorSet & solution);
private:
    Eigen::MatrixXd m_SubspaceMatrix;
};
}

#endif // DAVIDSON_H
