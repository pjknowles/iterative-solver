#include "RSPT.h"
using namespace IterativeSolver;

RSPT::RSPT(const ParameterSetTransformation residualFunction, const ParameterSetTransformation preconditionerFunction)
  : IterativeSolverBase(residualFunction, preconditionerFunction)
  , m_roots(1)
{
  m_linear = true;
  m_orthogonalize = false;
}

void RSPT::extrapolate(ParameterVectorSet & residual, ParameterVectorSet & solution, ParameterVectorSet & other, const optionMap options)
{
  size_t n=m_solutions.size();
  // on entry, solution contains |n-1> and residual contains H|n-1>, already stored in m_solutions & m_residuals.
  // on exit, incremental_energies[n] contains E_n = <0 |H-H_0-E_1|n-1> = <0|H|n-1>-(E_0+E_1)<0|n-1>.
  // eventually, Wigner 2n+1 rule should be implemented.
  // on exit, residual contains -(H_0-E_0)|n> = (H-H_0-E_1)|n-1> - sum_{k=0}^{n-2} E_{n-k} |k>.
  // = (H-E_0-E_1)|n-1> - (H_0-E_0)|n-1> - sum_{k=0}^{n-2} E_{n-k} |k>.
  // on exit, solution contains 0.

  diagonalizeSubspaceMatrix();
  if (m_verbosity>1) xout << "Subspace eigenvalues"<<std::endl<<m_subspaceEigenvalues<<std::endl;
  if (m_verbosity>2) xout << "Subspace eigenvectors"<<std::endl<<m_subspaceEigenvectors<<std::endl;

  if (n == 1) {
      std::vector<double> shift; shift.push_back(0);
      m_preconditionerFunction(solution,residual,shift,false);
      m_E0 = -1/solution.front()->dot(residual.front());
      solution.zero();
      m_incremental_energies.resize(2);
      m_incremental_energies[0]=m_E0;
      m_incremental_energies[1]=m_subspaceMatrix(0,0)-m_incremental_energies[0];
      residual.axpy(-m_subspaceMatrix(0,0),m_solutions.back());
      m_updateShift.clear();m_updateShift.push_back(-m_E0);
  } else {
      m_incremental_energies.push_back(m_subspaceMatrix(n-1,0)-m_subspaceMatrix(0,0)*m_subspaceOverlap(n-1,0));
      residual.axpy(1,m_lastH0mE0psi);
      residual.axpy(-m_subspaceMatrix(0,0),m_solutions.back());
      solution.zero();
  // this is structured for multistate, but not thought about yet
  for (size_t kkk=0; kkk<residual.size(); kkk++) {
      size_t l=0;
      for (size_t ll=0; ll<m_solutions.size(); ll++) {
          for (size_t lll=0; lll<m_solutions[ll].size(); lll++) {
//                  solution[kkk]->axpy(m_subspaceEigenvectors(l,kkk).real(),m_solutions[ll][lll]);
                  residual[kkk]->axpy(m_subspaceEigenvectors(l,kkk).real(),m_residuals[ll][lll]);
                  l++;
            }
        }
      residual[kkk]->axpy(-m_subspaceEigenvalues(kkk).real(),solution[kkk]);
    }

  }
  m_updateShift.resize(m_roots);
  for (size_t root=0; root<(size_t)m_roots; root++) m_updateShift[root]=m_singularity_shift-m_subspaceEigenvalues[root].real();
  m_lastH0mE0psi = residual; // we will need this in the next iteration // FIXME does this leak memory?
}
