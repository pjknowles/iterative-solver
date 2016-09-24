#include "IterativeSolver.h"
#include <algorithm>
#include <numeric>
#include <stdexcept>

using namespace IterativeSolver;

IterativeSolverBase::IterativeSolverBase(ParameterSetTransformation preconditionerFunction, ParameterSetTransformation residualFunction)
  :  m_preconditionerFunction(preconditionerFunction),
    m_residualFunction(residualFunction),
    m_verbosity(0),
    m_thresh(1e-12),
    m_maxIterations(1000),
    m_orthogonalize(false),
    m_linear(false),
    m_hermitian(false),
    m_preconditionResiduals(false),
    m_singularity_shift(1e-8)
{}

IterativeSolverBase::~IterativeSolverBase()
{
}


bool IterativeSolverBase::iterate(ParameterVectorSet &residual, ParameterVectorSet &solution, ParameterVectorSet &other, std::string options)
{
  if (m_preconditionResiduals) m_preconditionerFunction(residual,residual,m_updateShift,false);
  m_residuals.push_back(residual);
  m_solutions.push_back(solution); m_others.push_back(other);
  //  std::cout << "@@ solution at "<<&solution[0]<<std::endl;
  //  std::cout << "@@ m_solutions.back() at "<<&m_solutions.back()[0]<<std::endl;
  m_lastVectorIndex=m_residuals.size()-1; // derivative classes might eventually store the vectors on top of previous ones, in which case they will need to store the position here for later calculation of iteration step
  extrapolate(residual,solution,other,options);
  m_preconditionerFunction(residual,solution,m_updateShift,true);
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
  for (int iteration=1; iteration <= m_maxIterations && not converged; iteration++) {
      m_residualFunction(solution,residual,std::vector<double>(),false);
      converged = iterate(residual,solution);
      if (m_verbosity>0)
        xout << "iteration "<<iteration<<", error["<<m_worst<<"] = "<<m_error <<std::endl;
    }
  return converged;
}

void IterativeSolverBase::adjustUpdate(ParameterVectorSet &solution)
{
  for (size_t k=0; k<solution.size(); k++)
    solution.active[k] = (m_errors[k] > m_thresh);
  if (m_orthogonalize) {
      //      xout << "IterativeSolverBase::adjustUpdate solution before orthogonalization: "<<solution<<std::endl;
      for (size_t kkk=0; kkk<solution.size(); kkk++) {
          if (solution.active[kkk]) {
              size_t l=0;
              for (size_t ll=0; ll<m_solutions.size(); ll++) {
                  for (size_t lll=0; lll<m_solutions[ll].size(); lll++) {
                      if (m_solutions[ll].active[lll]) {
                          double s = -(m_solutions[ll][lll]->dot(solution[kkk])) / (m_solutions[ll][lll]->dot(m_solutions[ll][lll]));
                          solution[kkk]->axpy(s,m_solutions[ll][lll]);
                          l++;
                        }
                    }
                }
              for (size_t lll=0; lll<kkk; lll++) {
                  if (solution.active[lll]) {
                      double s = solution[lll]->dot(solution[kkk]);
                      solution[kkk]->axpy(-s,solution[lll]);
                    }
                }
              double s= solution[kkk]->dot(solution[kkk]);
              if (s <= 0)
                solution.active[kkk]=false;
              else
                solution[kkk]->axpy(1/std::sqrt(s)-1,solution[kkk]);
            }
        }
      //      xout << "IterativeSolverBase::adjustUpdate solution after orthogonalization: "<<solution<<std::endl;
    }
}

void IterativeSolverBase::calculateSubspaceMatrix(ParameterVectorSet &residual, ParameterVectorSet &solution)
{
  size_t old_size=m_subspaceMatrix.rows();
  m_subspaceMatrix.conservativeResize(old_size+residual.size(),old_size+residual.size());
  m_subspaceOverlap.conservativeResize(old_size+residual.size(),old_size+residual.size());
  size_t k=old_size;
  for (size_t kkk=0; kkk<residual.size(); kkk++) {
      if (residual.active[kkk]) {
          size_t l=0;
          for (size_t ll=0; ll<m_solutions.size(); ll++) {
              for (size_t lll=0; lll<m_solutions[ll].size(); lll++) {
                  if (m_solutions[ll].active[lll]) {
                      m_subspaceMatrix(k,l) = m_subspaceMatrix(l,k) = m_solutions[ll][lll]->dot(residual[kkk]);
                      m_subspaceOverlap(k,l) = m_subspaceOverlap(l,k) = m_solutions[ll][lll]->dot(solution[kkk]);
                      //                  std::cout << "subspace calc, *m_solutions[ll][lll] "<<*m_solutions[ll][lll]<<std::endl;
                      //                  std::cout << "subspace calc, solution[kkk] "<<*solution[kkk]<<std::endl;
                      //                  std::cout << "subspace calc, solution[kkk] "<<solution[kkk]->str()<<std::endl;
                      //                  std::cout << "subspace calc, product "<<solution[kkk]->dot(m_solutions[ll][lll])<<std::endl;
                      //                  std::cout << "subspace calc, product "<<m_solutions[ll][lll]->dot(solution[kkk])<<std::endl;
                      l++;
                    }
                }
            }
          k++;
        }
    }
  if (m_verbosity>3) {
      xout << "m_subspaceMatrix: "<<m_subspaceMatrix<<std::endl;
      xout << "m_subspaceOverlap: "<<m_subspaceOverlap<<std::endl;
    }

}

void IterativeSolverBase::diagonalizeSubspaceMatrix()
{
  Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> s(m_subspaceMatrix,m_subspaceOverlap);
  m_subspaceEigenvalues=s.eigenvalues();
  m_subspaceEigenvectors=s.eigenvectors();
  // sort
  std::vector<size_t> map;
  for (Eigen::Index k=0; k<m_subspaceMatrix.rows(); k++) {
      size_t ll;
      for (ll=0; std::count(map.begin(),map.end(),ll)!=0; ll++) ;
      for (Eigen::Index l=0; l<m_subspaceMatrix.rows(); l++) {
          if (std::count(map.begin(),map.end(),l)==0) {
              if (s.eigenvalues()(l).real() < s.eigenvalues()(ll).real())
                ll=l;
            }
        }
      map.push_back(ll);
      m_subspaceEigenvalues[k]=s.eigenvalues()(ll);
      for (Eigen::Index l=0; l<m_subspaceMatrix.rows(); l++) m_subspaceEigenvectors(l,k)=s.eigenvectors()(l,ll);
    }
  Eigen::MatrixXcd overlap=m_subspaceEigenvectors.transpose()*m_subspaceOverlap*m_subspaceEigenvectors;
  for (Eigen::Index k=0; k<overlap.rows(); k++)
    for (Eigen::Index l=0; l<overlap.rows(); l++)
      m_subspaceEigenvectors(l,k) /= std::sqrt(overlap(k,k).real());
}


void IterativeSolverBase::calculateErrors(const ParameterVectorSet &solution, const ParameterVectorSet &residual)
{
  ParameterVectorSet step=solution;
  step.axpy(-1,m_solutions[m_lastVectorIndex]);
  m_errors.clear();
  for (size_t k=0; k<solution.size(); k++)
    if (m_linear) // we can use the extrapolated residual if the problem is linear
      m_errors.push_back(residual.active[k] ? std::fabs(residual[k]->dot(step[k])) : 0);
    else
      m_errors.push_back(m_residuals[m_lastVectorIndex].active[k] ? std::fabs(m_residuals[m_lastVectorIndex][k]->dot(step[k])) : 0);
  m_error = *max_element(m_errors.begin(),m_errors.end());
  m_worst = max_element(m_errors.begin(),m_errors.end())-m_errors.begin();
}
