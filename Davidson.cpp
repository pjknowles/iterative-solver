#include "Davidson.h"
using namespace IterativeSolver;

Davidson::Davidson(ParameterSetTransformation updateFunction, ParameterSetTransformation residualFunction)
  : IterativeSolverBase(updateFunction, residualFunction)
  , m_roots(-1)
{
}

void Davidson::extrapolate(ParameterVectorSet & residual, ParameterVectorSet & solution, ParameterVectorSet & other, std::string options)
{
  assert(solution.size()==residual.size());
  std::cout << "residual.size() "<<residual.size()<<std::endl;
  std::cout << "solution.size() "<<solution.size()<<std::endl;
  if (m_roots<1) m_roots=solution.size(); // number of roots defaults to size of solution
  size_t old_size=m_SubspaceMatrix.rows();
  m_SubspaceMatrix.conservativeResize(old_size+residual.size(),old_size+residual.size());
  for (size_t k=0; k<residual.size(); k++) {
      size_t l=0;
      for (size_t ll=0; ll<m_solutions.size(); ll++) {
          for (size_t lll=0; lll<m_solutions[ll].size(); lll++) {
              if (m_solutions[ll].active[lll]) {
                  m_SubspaceMatrix(old_size+k,l) = m_SubspaceMatrix(l,old_size+k) = m_solutions[ll][lll] * residual[k];
                  m_SubspaceOverlap(old_size+k,l) = m_SubspaceOverlap(l,old_size+k) = m_solutions[ll][lll] * solution[k];
                  l++;
                }
            }
        }
    }
  std::cout << "Davidson::extrapolate m_SubspaceMatrix: "<<m_SubspaceMatrix<<std::endl;
  Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> s(m_SubspaceMatrix,m_SubspaceOverlap);
  m_Eigenvalues=s.eigenvalues();
  m_SubspaceEigenvectors=s.eigenvectors();
  std::cout << "eigenvalues"<<std::endl<<m_Eigenvalues<<std::endl;
  std::cout << "eigenvectors"<<std::endl<<m_SubspaceEigenvectors<<std::endl;
  m_updateShift.resize(m_roots);
  for (size_t root=0; root<m_roots; root++) m_updateShift[root]=1e-10-m_Eigenvalues[root];
}


static Eigen::MatrixXd testmatrix;

static void _residual(const ParameterVectorSet & psx, ParameterVectorSet & outputs, std::vector<ParameterScalar> shift=std::vector<ParameterScalar>()) {
  for (size_t k=0; k<psx.size(); k++) {
    Eigen::VectorXd x(testmatrix.rows());
    for (size_t l=0; l<testmatrix.rows(); l++) x[l] = psx[k][l];
    Eigen::VectorXd res = testmatrix * x;
    for (size_t l=0; l<testmatrix.rows(); l++) outputs[k][l]=res[l];
    }
}

static void _updater(const ParameterVectorSet & psg, ParameterVectorSet & psc, std::vector<ParameterScalar> shift) {
  std::cout << "updater: shift.size()="<<shift.size()<<std::endl;
  std::cout << "updater: psc="<<psc<<std::endl;
  std::cout << "updater: psg="<<psg<<std::endl;
  for (size_t k=0; k<psc.size(); k++) {
      std::cout << "shift "<<shift[k]<<std::endl;
    for (size_t l=0; l<testmatrix.rows(); l++)  psc[k][l] -= psg[k][l]/(testmatrix(l,l)+shift[k]);
    }
  std::cout << "updater: psc="<<psc<<std::endl;
}


void Davidson::test(size_t dimension, size_t roots, int verbosity)
{
  testmatrix.resize(dimension,dimension);
  for (size_t k=0; k<dimension; k++)
    for (size_t l=0; l<dimension; l++)
      testmatrix(l,k)=1;

  Davidson d(&_updater,&_residual);
  d.m_roots=roots;
  d.m_verbosity=verbosity;
  d.m_maxIterations=2;
  ParameterVectorSet x;
  ParameterVectorSet g;
  for (size_t root=0; root<d.m_roots; root++) {
      x.push_back(ParameterVector(dimension));
      x.front().zero();x.front()[root]=1;
      g.push_back(ParameterVector(dimension));
    }
  std::cout << "roots="<<roots<<std::endl;
  std::cout << "initial x="<<x<<std::endl;

  d.solve(g,x);


}
