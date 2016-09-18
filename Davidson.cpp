#include "Davidson.h"
using namespace IterativeSolver;

Davidson::Davidson(ParameterSetTransformation updateFunction, ParameterSetTransformation residualFunction)
  : IterativeSolverBase(updateFunction, residualFunction)
  , m_roots(-1)
{
  m_orthogonalize = true;
  m_true_extrapolated_residual = true;
}

void Davidson::extrapolate(ParameterVectorSet & residual, ParameterVectorSet & solution, ParameterVectorSet & other, std::string options)
{
  assert(solution.size()==residual.size());
  std::cout << "entry to extrapolate, residual "<<residual<<std::endl;
  std::cout << "entry to extrapolate, solution "<<solution<<std::endl;
  if (m_roots<1) m_roots=solution.size(); // number of roots defaults to size of solution
  size_t old_size=m_SubspaceMatrix.rows();
  m_SubspaceMatrix.conservativeResize(old_size+residual.size(),old_size+residual.size());
  m_SubspaceOverlap.conservativeResize(old_size+residual.size(),old_size+residual.size());
  size_t k=old_size;
  for (size_t kkk=0; kkk<residual.size(); kkk++) {
      if (residual.active[kkk]) {
      size_t l=0;
      for (size_t ll=0; ll<m_solutions.size(); ll++) {
          for (size_t lll=0; lll<m_solutions[ll].size(); lll++) {
              if (m_solutions[ll].active[lll]) {
                  m_SubspaceMatrix(k,l) = m_SubspaceMatrix(l,k) = m_solutions[ll][lll] * residual[kkk];
                  m_SubspaceOverlap(k,l) = m_SubspaceOverlap(l,k) = m_solutions[ll][lll] * solution[kkk];
                  l++;
                }
            }
        }
      k++;
      }
    }
  std::cout << "Davidson::extrapolate m_SubspaceMatrix: "<<m_SubspaceMatrix<<std::endl;
  std::cout << "Davidson::extrapolate m_SubspaceOverlap: "<<m_SubspaceOverlap<<std::endl;
  {
  Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> s(m_SubspaceMatrix,m_SubspaceOverlap);
  m_Eigenvalues=s.eigenvalues();
  std::cout << "Unsorted eigenvalues: "<<s.eigenvalues()<<std::endl;
  m_SubspaceEigenvectors=s.eigenvectors();
  // sort
  std::vector<size_t> map;
  for (size_t k=0; k<m_SubspaceMatrix.rows(); k++) {
      size_t ll;
      for (ll=0; std::count(map.begin(),map.end(),ll)!=0; ll++) ;
      for (size_t l=0; l<m_SubspaceMatrix.rows(); l++) {
          if (std::count(map.begin(),map.end(),l)==0) {
              if (s.eigenvalues()(l).real() < s.eigenvalues()(ll).real())
                  ll=l;
          }
      }
      map.push_back(ll);
      m_Eigenvalues[k]=s.eigenvalues()(ll);
      for (size_t l=0; l<m_SubspaceMatrix.rows(); l++) m_SubspaceEigenvectors(l,k)=s.eigenvectors()(l,ll);
  }
  Eigen::MatrixXcd overlap=m_SubspaceEigenvectors.transpose()*m_SubspaceOverlap*m_SubspaceEigenvectors;
  for (size_t k=0; k<overlap.rows(); k++)
      for (size_t l=0; l<overlap.rows(); l++)
          m_SubspaceEigenvectors(l,k) /= std::sqrt(overlap(k,k).real());
  }


  std::cout << "eigenvalues"<<std::endl<<m_Eigenvalues<<std::endl;
  std::cout << "eigenvectors"<<std::endl<<m_SubspaceEigenvectors<<std::endl;
  residual.zero();
  solution.zero();
  for (size_t kkk=0; kkk<residual.size(); kkk++) {
      size_t l=0;
      for (size_t ll=0; ll<m_solutions.size(); ll++) {
          for (size_t lll=0; lll<m_solutions[ll].size(); lll++) {
              if (m_solutions[ll].active[lll]) {
                  solution[kkk].axpy(m_SubspaceEigenvectors(l,kkk).real(),m_solutions[ll][lll]);
                  residual[kkk].axpy(m_SubspaceEigenvectors(l,kkk).real(),m_residuals[ll][lll]);
                  l++;
              }
          }
      }
      residual[kkk].axpy(-m_Eigenvalues(kkk).real(),solution[kkk]);
  }

  m_updateShift.resize(m_roots);
  for (size_t root=0; root<m_roots; root++) m_updateShift[root]=1e-10-m_Eigenvalues[root].real();

  std::cout << "exit from extrapolate, residual "<<residual<<std::endl;
  std::cout << "exit from extrapolate, solution "<<solution<<std::endl;
}

Eigen::VectorXd Davidson::eigenvalues()
{
Eigen::VectorXd result(m_roots);
  for (size_t root=0; root<m_roots; root++) result(root)=m_Eigenvalues[root].real();
  return result;
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
      testmatrix(l,k)=-1;

  Davidson d(&_updater,&_residual);
  d.m_roots=roots;
  d.m_verbosity=verbosity;
  d.m_maxIterations=dimension;
  ParameterVectorSet x;
  ParameterVectorSet g;
  for (size_t root=0; root<d.m_roots; root++) {
      x.push_back(ParameterVector(dimension));
      x.back().zero();x.back()[root]=1;
      std::cout << "x.size() "<<x.size()<<std::endl;
      std::cout << "x.back() "<<x.back()<<std::endl;
      g.push_back(ParameterVector(dimension));
    }
  std::cout << "roots="<<roots<<std::endl;
  std::cout << "initial x="<<x<<std::endl;

  d.solve(g,x);

  std::cout << "Eigenvalues: "<<d.eigenvalues()<<std::endl;


}
