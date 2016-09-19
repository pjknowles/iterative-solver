#include "IterativeSolver.h"
#include <algorithm>
#include <numeric>
#include <stdexcept>

using namespace IterativeSolver;

IterativeSolverBase::IterativeSolverBase(ParameterSetTransformation updateFunction, ParameterSetTransformation residualFunction)
:  m_updateFunction(updateFunction), m_residualFunction(residualFunction), m_verbosity(0), m_thresh(1e-12), m_maxIterations(1000), m_orthogonalize(false)
, m_linear(false), m_hermitian(false), m_singularity_shift(1e-8)
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
                  m_subspaceMatrix(k,l) = m_subspaceMatrix(l,k) = m_solutions[ll][lll] * residual[kkk];
                  m_subspaceOverlap(k,l) = m_subspaceOverlap(l,k) = m_solutions[ll][lll] * solution[kkk];
                  l++;
                }
            }
        }
      k++;
      }
    }
  if (m_verbosity>3) {
      std::cout << "m_subspaceMatrix: "<<m_subspaceMatrix<<std::endl;
      std::cout << "m_subspaceOverlap: "<<m_subspaceOverlap<<std::endl;
    }

}

void IterativeSolverBase::diagonalizeSubspaceMatrix()
{
  Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> s(m_subspaceMatrix,m_subspaceOverlap);
  m_subspaceEigenvalues=s.eigenvalues();
  m_subspaceEigenvectors=s.eigenvectors();
  // sort
  std::vector<size_t> map;
  for (size_t k=0; k<m_subspaceMatrix.rows(); k++) {
      size_t ll;
      for (ll=0; std::count(map.begin(),map.end(),ll)!=0; ll++) ;
      for (size_t l=0; l<m_subspaceMatrix.rows(); l++) {
          if (std::count(map.begin(),map.end(),l)==0) {
              if (s.eigenvalues()(l).real() < s.eigenvalues()(ll).real())
                  ll=l;
          }
      }
      map.push_back(ll);
      m_subspaceEigenvalues[k]=s.eigenvalues()(ll);
      for (size_t l=0; l<m_subspaceMatrix.rows(); l++) m_subspaceEigenvectors(l,k)=s.eigenvectors()(l,ll);
  }
  Eigen::MatrixXcd overlap=m_subspaceEigenvectors.transpose()*m_subspaceOverlap*m_subspaceEigenvectors;
  for (size_t k=0; k<overlap.rows(); k++)
      for (size_t l=0; l<overlap.rows(); l++)
          m_subspaceEigenvectors(l,k) /= std::sqrt(overlap(k,k).real());
  }


void IterativeSolverBase::calculateErrors(const ParameterVectorSet &solution, const ParameterVectorSet &residual)
{
  ParameterVectorSet step=solution;
  step.axpy(-1,m_solutions[m_lastVectorIndex]);
  m_errors.clear();
  for (size_t k=0; k<solution.size(); k++)
    if (m_linear) // we can use the extrapolated residual if the problem is linear
      m_errors.push_back(residual.active[k] ? std::fabs(residual[k]*step[k]) : 0);
    else
      m_errors.push_back(m_residuals[m_lastVectorIndex].active[k] ? std::fabs(m_residuals[m_lastVectorIndex][k]*step[k]) : 0);
  m_error = *max_element(m_errors.begin(),m_errors.end());
  m_worst = max_element(m_errors.begin(),m_errors.end())-m_errors.begin();
}

#if defined(__PGI) && ! defined(MOLPRO)
// https://www.molpro.net/bugzilla/show_bug.cgi?id=5012
__m128d _mm_castsi128_pd(__m128i x) { __m128d result; memcpy(&result,&x,16); }
__m128i _mm_castpd_si128(__m128d x) { __m128i result; memcpy(&result,&x,16); }
#endif
