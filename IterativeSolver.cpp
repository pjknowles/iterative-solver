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
    m_date(0),
    m_subspaceMatrixResRes(false),
    m_singularity_shift(1e-8)
{}

IterativeSolverBase::~IterativeSolverBase()
{
}


bool IterativeSolverBase::iterate(ParameterVectorSet &residual, ParameterVectorSet &solution, ParameterVectorSet &other, std::string options)
{
  for (size_t k=0; k<residual.size(); k++) residual.m_active[k] = residual.m_active[k] && solution.m_active[k];
  if (m_preconditionResiduals) m_preconditionerFunction(residual,residual,m_updateShift,false);
  m_lastVectorIndex=addVectorSet(residual,solution,other)-1; // derivative classes might eventually store the vectors on top of previous ones, in which case they will need to store the position here for later calculation of iteration step
  extrapolate(residual,solution,other,options);
  m_preconditionerFunction(residual,solution,m_updateShift,true);
  calculateErrors(solution,residual);
//  xout << "Errors:"; for (auto e=m_errors.begin(); e!=m_errors.end(); e++) xout << " "<<*e<<"("<<(*e<m_thresh)<<")"; xout <<std::endl;
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
    solution.m_active[k] = (m_errors[k] > m_thresh);
  if (m_orthogonalize) {
      //      xout << "IterativeSolverBase::adjustUpdate solution before orthogonalization: "<<solution<<std::endl;
      for (size_t kkk=0; kkk<solution.size(); kkk++) {
          if (solution.m_active[kkk]) {
              size_t l=0;
              for (size_t ll=0; ll<m_solutions.size(); ll++) {
                  for (size_t lll=0; lll<m_solutions[ll].size(); lll++) {
                      if (m_solutions[ll].m_active[lll]) {
                          double s = -(m_solutions[ll][lll]->dot(solution[kkk])) / (m_solutions[ll][lll]->dot(m_solutions[ll][lll]));
                          solution[kkk]->axpy(s,m_solutions[ll][lll]);
                          l++;
                        }
                    }
                }
              for (size_t lll=0; lll<kkk; lll++) {
                  if (solution.m_active[lll]) {
                      double s = solution[lll]->dot(solution[kkk]);
                      solution[kkk]->axpy(-s,solution[lll]);
                    }
                }
              double s= solution[kkk]->dot(solution[kkk]);
              if (s <= 0)
                solution.m_active[kkk]=false;
              else
                solution[kkk]->axpy(1/std::sqrt(s)-1,solution[kkk]);
            }
        }
      //      xout << "IterativeSolverBase::adjustUpdate solution after orthogonalization: "<<solution<<std::endl;
    }
}

size_t IterativeSolverBase::addVectorSet(const ParameterVectorSet &residual, const ParameterVectorSet &solution, const ParameterVectorSet &other)
{
    m_residuals.push_back(residual);
    m_solutions.push_back(solution);
    m_others.push_back(other);
    calculateSubspaceMatrix(residual,solution);
    return m_residuals.size();
}

void IterativeSolverBase::deleteVector(size_t index)
{
    if (index>=m_residuals.size()) throw std::logic_error("invalid index");
    xout << "deleteVector "<<index<<std::endl;
    xout << "old m_subspaceMatrix"<<std::endl<<m_subspaceMatrix<<std::endl;
    xout << "old m_subspaceOverlap"<<std::endl<<m_subspaceOverlap<<std::endl;
    size_t old_size=m_subspaceMatrix.rows();
    size_t new_size=old_size;
    size_t l=0;
    for (size_t ll=0; ll<m_solutions.size(); ll++) {
        for (size_t lll=0; lll<m_solutions[ll].size(); lll++) {
            if (m_solutions[ll].m_active[lll]) {
                if (l == index) { // this is the one to delete
                    m_solutions[ll].m_active[lll] = false;
                    m_residuals[ll].m_active[lll] = false;
                    m_others[ll].m_active[lll] = false;
                    for (size_t l2=l+1; l2<m_subspaceMatrix.rows(); l2++) {
                        for (size_t k=0; k<m_subspaceMatrix.rows(); k++) {
                            m_subspaceMatrix(k,l2-1)=m_subspaceMatrix(k,l2);
                            m_subspaceMatrix(l2-1,k)=m_subspaceMatrix(l2,k);
                            m_subspaceOverlap(k,l2-1)=m_subspaceOverlap(k,l2);
                            m_subspaceOverlap(l2-1,k)=m_subspaceOverlap(l2,k);
                        }
                    }
                    new_size--;
                }
                l++;
            }
            if (std::count(m_solutions[ll].m_active.begin(),m_solutions[ll].m_active.end(),true)==0) { // can delete the whole thing
            }
        }
    }
    m_subspaceMatrix.conservativeResize(new_size,new_size);
    m_subspaceOverlap.conservativeResize(new_size,new_size);
    xout << "new m_subspaceMatrix"<<std::endl<<m_subspaceMatrix<<std::endl;
    xout << "new m_subspaceOverlap"<<std::endl<<m_subspaceOverlap<<std::endl;
}

void IterativeSolverBase::calculateSubspaceMatrix(const ParameterVectorSet &residual, const ParameterVectorSet &solution)
{
  size_t old_size=m_subspaceMatrix.rows();
  size_t new_size=old_size+std::count(residual.m_active.begin(),residual.m_active.end(),true);
  m_subspaceMatrix.conservativeResize(new_size,new_size);
  m_subspaceOverlap.conservativeResize(new_size,new_size);
  std::vector<ParameterVectorSet>* bra = m_subspaceMatrixResRes ? &m_residuals : &m_solutions;
  size_t k=old_size;
  for (size_t kkk=0; kkk<residual.size(); kkk++) {
      if (residual.m_active[kkk]) {
          size_t l=0;
          for (size_t ll=0; ll<m_solutions.size(); ll++) {
              for (size_t lll=0; lll<m_solutions[ll].size(); lll++) {
                  if (m_solutions[ll].m_active[lll]) {
                      m_subspaceMatrix(k,l) = m_subspaceMatrix(l,k) = (*bra)[ll][lll]->dot(residual[kkk]);
                      m_subspaceOverlap(k,l) = m_subspaceOverlap(l,k) = m_solutions[ll][lll]->dot(solution[kkk]);
                      l++;
                    }
                }
            }
          m_dateOfBirth.push_back(++m_date);
          k++;
        }
    }
  if (m_verbosity>3) {
      xout << "m_subspaceMatrix: "<<std::endl<<m_subspaceMatrix<<std::endl;
      xout << "m_subspaceOverlap: "<<std::endl<<m_subspaceOverlap<<std::endl;
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
      m_errors.push_back(residual.m_active[k] ? std::fabs(residual[k]->dot(step[k])) : 0);
    else
      m_errors.push_back(m_residuals[m_lastVectorIndex].m_active[k] ? std::fabs(m_residuals[m_lastVectorIndex][k]->dot(step[k])) : 0);
  m_error = *max_element(m_errors.begin(),m_errors.end());
  m_worst = max_element(m_errors.begin(),m_errors.end())-m_errors.begin();
}
