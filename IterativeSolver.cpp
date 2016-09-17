#include "IterativeSolver.h"
#include <algorithm>
#include <numeric>

using namespace IterativeSolver;

IterativeSolverBase::IterativeSolverBase(ParameterSetTransformation updateFunction, ParameterSetTransformation residualFunction)
:  m_updateFunction(updateFunction), m_residualFunction(residualFunction), m_verbosity(0), m_thresh(1e-12), m_maxIterations(1000)
{}

IterativeSolverBase::~IterativeSolverBase()
{
}


#include <limits>

bool IterativeSolverBase::iterate(ParameterVectorSet &residual, ParameterVectorSet &solution, ParameterVectorSet &other, std::string options)
{
	m_residuals.push_back(residual); m_solutions.push_back(solution); m_others.push_back(other);
	m_lastVectorIndex=m_residuals.size()-1; // derivative classes might eventually store the vectors on top of previous ones, in which case they will need to store the position here for later calculation of iteration step
    if (m_verbosity>3) {
	  std::cout <<"added to m_residuals, now size="<<m_residuals.size()<<std::endl;
      std::cout << "latest residual: "<<m_residuals.back()<<std::endl;
      std::cout << "latest solution: "<<m_solutions.back()<<std::endl;
    }
    extrapolate(residual,solution,other,options);
//	  std::cout << "after extrapolate solution: "<<solution<<std::endl;
	m_updateFunction(residual,solution,std::vector<ParameterScalar>());
    calculateErrors(solution);
//    std::cout << "iterate() final extrapolated and updated solution: "<<solution<<std::endl;
//    std::cout << "iterate() final extrapolated and updated residual: "<<residual<<std::endl;
//    std::cout << "error: "<<err<<", m_thresh="<<m_thresh<<std::endl;
        return m_error < m_thresh;
}


 void IterativeSolverBase::extrapolate(ParameterVectorSet & residual, ParameterVectorSet & solution, ParameterVectorSet & other, std::string options)
 {

 }

 bool IterativeSolverBase::solve(ParameterVectorSet & residual, ParameterVectorSet & solution, std::string options)
 {
         bool converged=false;
         size_t max_iterations=1000; //FIXME
     for (size_t iteration=1; iteration <= m_maxIterations && not converged; iteration++) {
         m_residualFunction(solution,residual,std::vector<double>());
         converged = iterate(residual,solution);
         if (m_verbosity>0)
         std::cout << "iteration "<<iteration<<", error["<<m_worst<<"] = "<<m_error
         <<std::endl;
     }
   return converged;
 }

void IterativeSolverBase::calculateErrors(const ParameterVectorSet &solution)
{
  ParameterVectorSet step=solution;
  step.axpy(-1,m_solutions[m_lastVectorIndex]);
//  std::cout << "m_solutions[m_lastVectorIndex]="<<m_solutions[m_lastVectorIndex]<<std::endl;
//  std::cout << "solution="<<solution<<std::endl;
//  std::cout << "m_residuals[m_lastVectorIndex]="<<m_residuals[m_lastVectorIndex]<<std::endl;
//  std::cout << "step="<<step<<std::endl;
  m_errors.clear();
    for (size_t k=0; k<solution.size(); k++)
        m_errors.push_back(m_residuals[m_lastVectorIndex].active[k] ? std::fabs(m_residuals[m_lastVectorIndex][k]*step[k]) : 0);
    m_error = *max_element(m_errors.begin(),m_errors.end());
    m_worst = max_element(m_errors.begin(),m_errors.end())-m_errors.begin();
}
