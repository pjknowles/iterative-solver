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
  assert(solution.size()==residual.size());
  if (m_roots<1) m_roots=solution.size(); // number of roots defaults to size of solution
//  calculateSubspaceMatrix(residual,solution);
  diagonalizeSubspaceMatrix();

  if (m_verbosity>1) xout << "Subspace eigenvalues"<<std::endl<<m_subspaceEigenvalues<<std::endl;
  if (m_verbosity>2) xout << "Subspace eigenvectors"<<std::endl<<m_subspaceEigenvectors<<std::endl;
  residual.zero();
  solution.zero();
  for (size_t kkk=0; kkk<residual.size(); kkk++) {
      size_t l=0;
      for (size_t ll=0; ll<m_solutions.size(); ll++) {
          for (size_t lll=0; lll<m_solutions[ll].size(); lll++) {
              if (m_solutions[ll].m_active[lll]) {
                  solution[kkk]->axpy(m_subspaceEigenvectors(l,kkk).real(),m_solutions[ll][lll]);
                  residual[kkk]->axpy(m_subspaceEigenvectors(l,kkk).real(),m_residuals[ll][lll]);
                  l++;
                }
            }
        }
      residual[kkk]->axpy(-m_subspaceEigenvalues(kkk).real(),solution[kkk]);
    }

  m_updateShift.resize(m_roots);
  for (size_t root=0; root<(size_t)m_roots; root++) m_updateShift[root]=m_singularity_shift-m_subspaceEigenvalues[root].real();
}
