#include "IterativeSolver.h"
#include <algorithm>
#include <numeric>
#include <stdexcept>

using namespace IterativeSolver;

IterativeSolverBase::IterativeSolverBase(const ParameterSetTransformation residualFunction, const ParameterSetTransformation preconditionerFunction)
  :  m_preconditionerFunction(preconditionerFunction),
    m_residualFunction(residualFunction),
    m_verbosity(0),
    m_thresh(1e-12),
    m_maxIterations(1000),
    m_minIterations(0),
    m_orthogonalize(false),
    m_linear(false),
    m_hermitian(false),
    m_preconditionResiduals(false),
    m_roots(-1),
    m_date(0),
    m_subspaceMatrixResRes(false),
    m_singularity_shift(1e-8),
    m_iterations(0),
    m_singularity_threshold(1e-20)
{}

IterativeSolverBase::~IterativeSolverBase()
{
}


bool IterativeSolverBase::iterate(ParameterVectorSet &residual, ParameterVectorSet &solution, ParameterVectorSet &other, const optionMap options)
{
//  xout << "iterate"<<std::endl;
//  xout << "residual"<<std::endl<<residual[0]<<std::endl;
//  xout << "solution"<<std::endl<<solution[0]<<std::endl;
  if (m_roots<1) m_roots=solution.size(); // number of roots defaults to size of solution
  assert(solution.size()==residual.size());
  m_iterations++;
  for (size_t k=0; k<residual.size(); k++) residual.m_active[k] = residual.m_active[k] && solution.m_active[k];
  if (m_preconditionResiduals) m_preconditionerFunction(residual,residual,m_updateShift,false);
  m_lastVectorIndex=addVectorSet(residual,solution,other)-1; // derivative classes might eventually store the vectors on top of previous ones, in which case they will need to store the position here for later calculation of iteration step
  extrapolate(residual,solution,other,options);
  if (m_preconditionResiduals)
    solution.axpy(1,residual);
  else
    m_preconditionerFunction(residual,solution,m_updateShift,true);
  calculateErrors(solution,residual);
//  xout << "Errors:"; for (auto e=m_errors.begin(); e!=m_errors.end(); e++) xout << " "<<*e<<"("<<(*e<m_thresh)<<")"; xout <<std::endl;
  adjustUpdate(solution);
  return m_error < m_thresh;
}


bool IterativeSolverBase::solve(ParameterVectorSet & residual, ParameterVectorSet & solution, const optionMap options)
{
  bool converged=false;
  for (unsigned int iteration=1; iteration <= m_maxIterations && (not converged || iteration <= m_minIterations); iteration++) {
      m_residualFunction(solution,residual,std::vector<double>(),false);
      converged = iterate(residual,solution);
      report();
    }
  return converged;
}

void IterativeSolverBase::report()
{
      if (m_verbosity>0)
        xout << "iteration "<<iterations()<<", error["<<m_worst<<"] = "<<m_error <<std::endl;
}

void IterativeSolverBase::adjustUpdate(ParameterVectorSet &solution)
{
  for (size_t k=0; k<solution.size(); k++)
    solution.m_active[k] = (m_errors[k] >= m_thresh || m_minIterations>m_iterations);
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
//  if (residual.m_active.front()==0) xout <<"warning: inactive residual"<<std::endl;
//  if (solution.m_active.front()==0) xout <<"warning: inactive solution"<<std::endl;
    m_residuals.push_back(residual);
    m_solutions.push_back(solution);
    m_others.push_back(other);
    calculateSubspaceMatrix(residual,solution);
    return m_residuals.size();
}

void IterativeSolverBase::deleteVector(size_t index)
{
    if (index>=m_residuals.size()) throw std::logic_error("invalid index");
//    xout << "deleteVector "<<index<<std::endl;
//    xout << "old m_subspaceMatrix"<<std::endl<<m_subspaceMatrix<<std::endl;
//    xout << "old m_subspaceOverlap"<<std::endl<<m_subspaceOverlap<<std::endl;
    size_t old_size=m_subspaceMatrix.rows();
    size_t new_size=old_size;
    size_t l=0;
    for (size_t ll=0; ll<m_solutions.size(); ll++) {
        for (size_t lll=0; lll<m_solutions[ll].size(); lll++) {
            if (m_solutions[ll].m_active[lll]) {
                if (l == index) { // this is the one to delete
                    m_dateOfBirth.erase(m_dateOfBirth.begin()+l);
                    m_solutions[ll].m_active[lll] = false;
                    m_residuals[ll].m_active[lll] = false;
//                    m_others[ll].m_active[lll] = false;
                    for (int l2=l+1; l2<m_subspaceMatrix.rows(); l2++) {
                        for (int k=0; k<m_subspaceMatrix.rows(); k++) {
//                            xout << "copy from ("<<k<<","<<l2<<") to ("<<k<<","<<l2-1<<")"<<std::endl;
                            m_subspaceMatrix(k,l2-1)=m_subspaceMatrix(k,l2);
                            m_subspaceOverlap(k,l2-1)=m_subspaceOverlap(k,l2);
                        }
                    }
                    for (int l2=l+1; l2<m_subspaceMatrix.rows(); l2++) {
                        for (int k=0; k<m_subspaceMatrix.rows(); k++) {
                            m_subspaceMatrix(l2-1,k)=m_subspaceMatrix(l2,k);
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
//    xout << "new m_subspaceMatrix"<<std::endl<<m_subspaceMatrix<<std::endl;
//    xout << "new m_subspaceOverlap"<<std::endl<<m_subspaceOverlap<<std::endl;
}

void IterativeSolverBase::calculateSubspaceMatrix(const ParameterVectorSet &residual, const ParameterVectorSet &solution)
{
//  xout << "calculateSubspaceMatrix"<<std::endl;
//  xout << "residual"<<std::endl<<residual[0]<<std::endl;
//  xout << "solution"<<std::endl<<solution[0]<<std::endl;
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
//  xout << "bra"<<std::endl<<(*bra)[ll][lll]<<std::endl;
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
  int kept=m_subspaceMatrix.rows();
  {
      Eigen::EigenSolver<Eigen::MatrixXd> ss(m_subspaceOverlap);
      Eigen::VectorXcd sse=ss.eigenvalues();
      for (int k=0; k<sse.rows(); k++) {
          if (std::fabs(sse(k).real()) < m_singularity_threshold)
              kept--;
      }
  }
  if (m_verbosity >=0 && kept < m_subspaceMatrix.rows())
      xout <<"IterativeSolver WARNING, subspace singular, pruned from "<<m_subspaceMatrix.rows()<<" to "<<kept<<std::endl;

  Eigen::MatrixXd H=m_subspaceMatrix.block(0,0,kept,kept);
  Eigen::MatrixXd S=m_subspaceOverlap.block(0,0,kept,kept);
  Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> s(H,S);
  m_subspaceEigenvalues=s.eigenvalues();
  m_subspaceEigenvectors=s.eigenvectors();
  // sort
  std::vector<size_t> map;
  for (Eigen::Index k=0; k<H.rows(); k++) {
      size_t ll;
      for (ll=0; std::count(map.begin(),map.end(),ll)!=0; ll++) ;
      for (Eigen::Index l=0; l<H.rows(); l++) {
          if (std::count(map.begin(),map.end(),l)==0) {
              if (s.eigenvalues()(l).real() < s.eigenvalues()(ll).real())
                ll=l;
            }
        }
      map.push_back(ll);
      m_subspaceEigenvalues[k]=s.eigenvalues()(ll);
      for (Eigen::Index l=0; l<H.rows(); l++) m_subspaceEigenvectors(l,k)=s.eigenvectors()(l,ll);
    }
  Eigen::MatrixXcd overlap=m_subspaceEigenvectors.transpose()*S*m_subspaceEigenvectors;
  for (Eigen::Index k=0; k<overlap.rows(); k++)
    for (Eigen::Index l=0; l<overlap.rows(); l++)
      m_subspaceEigenvectors(l,k) /= std::sqrt(overlap(k,k).real());
//  xout << "eigenvalues"<<std::endl<<m_subspaceEigenvalues<<std::endl;
//  xout << "eigenvectors"<<std::endl<<m_subspaceEigenvectors<<std::endl;
}


#include <math.h>
void IterativeSolverBase::calculateErrors(const ParameterVectorSet &solution, const ParameterVectorSet &residual)
{
  ParameterVectorSet step=solution;
  step.axpy(-1,m_solutions[m_lastVectorIndex]);
  m_errors.clear();
//  xout << "last active "<<m_lastVectorIndex<<" "<<m_residuals[m_lastVectorIndex].m_active[0]<<std::endl;
  for (size_t k=0; k<solution.size(); k++) {
    if (m_linear) // we can use the extrapolated residual if the problem is linear
      m_errors.push_back(residual.m_active[k] ? std::fabs(residual[k]->dot(step[k])) : 0);
    else
      m_errors.push_back(m_residuals[m_lastVectorIndex].m_active[k] ? std::fabs(m_residuals[m_lastVectorIndex][k]->dot(step[k])) : 1);
    if (isnan(m_errors.back())) throw std::overflow_error("NaN detected in DIIS error measure");
    }
  m_error = *max_element(m_errors.begin(),m_errors.end());
  m_worst = max_element(m_errors.begin(),m_errors.end())-m_errors.begin();
}

std::vector<double> IterativeSolverBase::eigenvalues()
{
  std::vector<double> result;
  for (unsigned int root=0; root<(size_t)m_roots && root < m_subspaceEigenvalues.rows(); root++) result.push_back(m_subspaceEigenvalues[root].real());
  return result;
}
