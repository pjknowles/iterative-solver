#include "IterativeSolver.h"
#include <algorithm>
#include <numeric>

using namespace IterativeSolver;

IterativeSolverBase::IterativeSolverBase(ParameterSetTransformation updateFunction, ParameterSetTransformation residualFunction)
:  m_updateFunction(updateFunction), m_residualFunction(residualFunction), m_verbosity(0), m_thresh(1e-12), m_maxIterations(1000), m_orthogonalize(false), m_true_extrapolated_residual(false)
{}

IterativeSolverBase::~IterativeSolverBase()
{
}


bool IterativeSolverBase::iterate(ParameterVectorSet &residual, ParameterVectorSet &solution, ParameterVectorSet &other, std::string options)
{
  m_residuals.push_back(residual); m_solutions.push_back(solution); m_others.push_back(other);
  m_lastVectorIndex=m_residuals.size()-1; // derivative classes might eventually store the vectors on top of previous ones, in which case they will need to store the position here for later calculation of iteration step
  extrapolate(residual,solution,other,options);
  m_updateFunction(residual,solution,m_updateShift);
  calculateErrors(solution,residual);
  adjustUpdate(solution);
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
        std::cout << "iteration "<<iteration<<", error["<<m_worst<<"] = "<<m_error <<std::endl;
    }
  return converged;
}

void IterativeSolverBase::adjustUpdate(ParameterVectorSet &solution)
{
  for (size_t k=0; k<solution.size(); k++)
    solution.active[k] = (m_errors[k] > m_thresh);
  if (m_orthogonalize) {
//      std::cout << "IterativeSolverBase::adjustUpdate solution before orthogonalization: "<<solution<<std::endl;
      for (size_t kkk=0; kkk<solution.size(); kkk++) {
          if (solution.active[kkk]) {
              size_t l=0;
              for (size_t ll=0; ll<m_solutions.size(); ll++) {
                  for (size_t lll=0; lll<m_solutions[ll].size(); lll++) {
                      if (m_solutions[ll].active[lll]) {
                          double s = -(m_solutions[ll][lll] * solution[kkk]) / (m_solutions[ll][lll] * m_solutions[ll][lll]);
                          solution[kkk].axpy(s,m_solutions[ll][lll]);
                          l++;
                        }
                    }
                }
              for (size_t lll=0; lll<kkk; lll++) {
                  if (solution.active[lll]) {
                      double s = solution[lll] * solution[kkk];
                      solution[kkk].axpy(-s,solution[lll]);
                    }
                }
              double s= solution[kkk]*solution[kkk];
              if (s <= 0)
                solution.active[kkk]=false;
              else
                solution[kkk].axpy(1/std::sqrt(s)-1,solution[kkk]);
            }
        }
//      std::cout << "IterativeSolverBase::adjustUpdate solution after orthogonalization: "<<solution<<std::endl;
    }
}

void IterativeSolverBase::calculateErrors(const ParameterVectorSet &solution, const ParameterVectorSet &residual)
{
  ParameterVectorSet step=solution;
  step.axpy(-1,m_solutions[m_lastVectorIndex]);
  m_errors.clear();
  for (size_t k=0; k<solution.size(); k++)
    if (m_true_extrapolated_residual)
      m_errors.push_back(residual.active[k] ? std::fabs(residual[k]*step[k]) : 0);
    else
      m_errors.push_back(m_residuals[m_lastVectorIndex].active[k] ? std::fabs(m_residuals[m_lastVectorIndex][k]*step[k]) : 0);
  m_error = *max_element(m_errors.begin(),m_errors.end());
  m_worst = max_element(m_errors.begin(),m_errors.end())-m_errors.begin();
}
