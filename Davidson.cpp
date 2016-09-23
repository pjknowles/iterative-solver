#include "Davidson.h"
#include <stdexcept>
using namespace IterativeSolver;

Davidson::Davidson(ParameterSetTransformation updateFunction, ParameterSetTransformation residualFunction)
  : IterativeSolverBase(updateFunction, residualFunction)
  , m_roots(-1)
{
  m_linear = true;
  m_orthogonalize = true;
}


void Davidson::extrapolate(ParameterVectorSet & residual, ParameterVectorSet & solution, ParameterVectorSet & other, std::string options)
{
  assert(solution.size()==residual.size());
//  xout << "entry to extrapolate, residual "<<residual<<std::endl;
//  xout << "entry to extrapolate, solution "<<solution<<std::endl;
  if (m_roots<1) m_roots=solution.size(); // number of roots defaults to size of solution
  calculateSubspaceMatrix(residual,solution);
  diagonalizeSubspaceMatrix();

  if (m_verbosity>1) xout << "Subspace eigenvalues"<<std::endl<<m_subspaceEigenvalues<<std::endl;
  if (m_verbosity>2) xout << "Subspace eigenvectors"<<std::endl<<m_subspaceEigenvectors<<std::endl;
  residual.zero();
  solution.zero();
  for (size_t kkk=0; kkk<residual.size(); kkk++) {
      size_t l=0;
      for (size_t ll=0; ll<m_solutions.size(); ll++) {
          for (size_t lll=0; lll<m_solutions[ll].size(); lll++) {
              if (m_solutions[ll].active[lll]) {
                  solution[kkk].axpy(m_subspaceEigenvectors(l,kkk).real(),m_solutions[ll][lll]);
                  residual[kkk].axpy(m_subspaceEigenvectors(l,kkk).real(),m_residuals[ll][lll]);
                  l++;
              }
          }
      }
      residual[kkk].axpy(-m_subspaceEigenvalues(kkk).real(),solution[kkk]);
  }

  m_updateShift.resize(m_roots);
  for (size_t root=0; root<(size_t)m_roots; root++) m_updateShift[root]=m_singularity_shift-m_subspaceEigenvalues[root].real();

//  xout << "exit from extrapolate, residual "<<residual<<std::endl;
//  xout << "exit from extrapolate, solution "<<solution<<std::endl;
}

std::vector<double> Davidson::eigenvalues()
{
std::vector<double> result;
  for (size_t root=0; root<(size_t)m_roots; root++) result.push_back(m_subspaceEigenvalues[root].real());
  return result;
}

// testing code below here
#include "SimpleParameterVector.h"
static Eigen::MatrixXd testmatrix;

static void _residual(const ParameterVectorSet & psx, ParameterVectorSet & outputs, std::vector<ParameterScalar> shift=std::vector<ParameterScalar>(), bool append=false) {
  for (size_t k=0; k<psx.size(); k++) {
    Eigen::VectorXd x(testmatrix.rows());
    for (size_t l=0; l<(size_t)testmatrix.rows(); l++) x[l] = psx[k][l];
    Eigen::VectorXd res = testmatrix * x;
    if (not append) outputs[k].zero();
    for (size_t l=0; l<(size_t)testmatrix.rows(); l++) outputs[k][l]+=res[l];
    }
}

static void _updater(const ParameterVectorSet & psg, ParameterVectorSet & psc, std::vector<ParameterScalar> shift, bool append=true) {
//  xout << "updater: shift.size()="<<shift.size()<<std::endl;
//  xout << "updater: psc="<<psc<<std::endl;
//  xout << "updater: psg="<<psg<<std::endl;
    if (not append) psc.zero();
  for (size_t k=0; k<psc.size(); k++) {
//      xout << "shift "<<shift[k]<<std::endl;
    for (size_t l=0; l<(size_t)testmatrix.rows(); l++)  psc[k][l] -= psg[k][l]/(testmatrix(l,l)+shift[k]);
    }
//  xout << "updater: psc="<<psc<<std::endl;
}


//typedef SimpleParameterVector ptype;
typedef ParameterVector ptype;
void Davidson::test(size_t dimension, size_t roots, int verbosity, int problem)
{
  xout << "Test IterativeSolver::Davidson dimension="<<dimension<<", roots="<<roots<<", problem="<<problem<<std::endl;
  testmatrix.resize(dimension,dimension);
  for (size_t k=0; k<dimension; k++)
    for (size_t l=0; l<dimension; l++)
      if (problem==0)
        testmatrix(l,k)=-1;
      else if (problem==1)
        testmatrix(l,k)=l+k+2;
      else if (problem==2)
        testmatrix(l,k)=( k==l ? k+1 : 1);
      else
        throw std::logic_error("invalid problem in Davidson::test");

  Davidson d(&_updater,&_residual);
  d.m_roots=roots;
  d.m_verbosity=verbosity;
  d.m_maxIterations=dimension;
  ParameterVectorSet x;
  ParameterVectorSet g;
  for (size_t root=0; root<(size_t)d.m_roots; root++) {
      ptype* xx=new ptype(dimension);
      xx->zero();
      (*xx)[root]=1;
      x.push_back(*xx);
      ptype* gg=new ptype(dimension);
      g.push_back(*gg);
    }
  xout << "roots="<<roots<<std::endl;
  xout << "initial x="<<x<<std::endl;

  d.solve(g,x);

  std::vector<double> ev=d.eigenvalues();
  xout << "Eigenvalues: "; for (std::vector<double>::const_iterator e=ev.begin(); e!=ev.end(); e++) xout<<" "<<*e;xout<<std::endl;
  xout << "Reported errors: "; for (std::vector<double>::const_iterator e=d.m_errors.begin(); e!=d.m_errors.end(); e++) xout<<" "<<*e;xout<<std::endl;

  _residual(g,x);
  std::vector<double> errors;
  for (size_t root=0; root<(size_t)d.m_roots; root++) {
      g[root].axpy(-ev[root],x[root]);
      errors.push_back(g[root]*g[root]);
    }
  xout << "Square residual norms: "; for (std::vector<double>::const_iterator e=errors.begin(); e!=errors.end(); e++) xout<<" "<<*e;xout<<std::endl;
  // be noisy about obvious problems
  if (*std::max_element(errors.begin(),errors.end())>1e-8) throw std::runtime_error("IterativeSolver::Davidson has failed tests");

  for (size_t root=0; root<(size_t)d.m_roots; root++) {
      delete &x[root][0];
      delete &g[root][0];
  }

}
