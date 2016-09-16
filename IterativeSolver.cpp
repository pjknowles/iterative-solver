#include "IterativeSolver.h"
#include <algorithm>
#include <numeric>

using namespace IterativeSolver;

IterativeSolverBase::IterativeSolverBase(ParameterSetTransformation updateFunction, ParameterSetTransformation residualFunction)
:  m_updateFunction(updateFunction), m_residualFunction(residualFunction), m_verbosity(0), m_thresh(1e-6)
{}

IterativeSolverBase::~IterativeSolverBase()
{
}


#include <limits>

bool IterativeSolverBase::iterate(ParameterVectorSet &residual, ParameterVectorSet &solution, ParameterVectorSet &other, std::string options)
{
	m_residuals.push_back(residual); m_solutions.push_back(solution); m_others.push_back(other);
    if (m_verbosity>3) {
	  std::cout <<"added to m_residuals, now size="<<m_residuals.size()<<std::endl;
      std::cout << "latest residual: "<<m_residuals.back()<<std::endl;
      std::cout << "latest solution: "<<m_solutions.back()<<std::endl;
    }
       double err = calculateError(residual,solution);
    extrapolate(residual,solution,other,options);
//	  std::cout << "after extrapolate solution: "<<solution<<std::endl;
	m_updateFunction(residual,solution,std::vector<ParameterScalar>());
//    std::cout << "iterate() final extrapolated and updated solution: "<<solution<<std::endl;
//    std::cout << "iterate() final extrapolated and updated residual: "<<residual<<std::endl;
//    std::cout << "error: "<<err<<", m_thresh="<<m_thresh<<std::endl;
        return err < m_thresh;
}


 void IterativeSolverBase::extrapolate(ParameterVectorSet & residual, ParameterVectorSet & solution, ParameterVectorSet & other, std::string options)
 {

 }

 bool IterativeSolverBase::solve(ParameterVectorSet & residual, ParameterVectorSet & solution, std::string options)
 {
	 throw std::logic_error("solve not written yet"); // FIXME
   return false;
 }

std::vector<double> IterativeSolverBase::calculateErrors(const ParameterVectorSet &residual, const ParameterVectorSet &solution)
{
    std::vector<double> result;
    for (size_t k=0; k<residual.size(); k++)
        result.push_back(residual.active[k] ? residual[k]*solution[k] : 0);
    return result;
}

double IterativeSolverBase::calculateError(const ParameterVectorSet &residual, const ParameterVectorSet &solution)
{
    std::vector<double> errs;
    errs=calculateErrors(residual,solution);
    double result=0;
    for (std::vector<double>::iterator e=errs.begin(); e!=errs.end(); e++) result+=(*e)*(*e);
    return std::sqrt(result);
}
