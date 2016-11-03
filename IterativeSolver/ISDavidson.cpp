#include "ISDavidson.h"
#include <stdexcept>
using namespace LinearAlgebra;

Davidson::Davidson(const ParameterSetTransformation residualFunction, const ParameterSetTransformation preconditionerFunction)
  : IterativeSolverBase(residualFunction, preconditionerFunction)
{
  m_linear = true;
  m_orthogonalize = true;
}


void Davidson::extrapolate(ParameterVectorSet & residual, ParameterVectorSet & solution, ParameterVectorSet & other, const optionMap options)
{
  if (m_verbosity>2) xout << "Subspace matrix"<<std::endl<<m_subspaceMatrix<<std::endl;
  if (m_verbosity>2) xout << "Subspace overlap"<<std::endl<<m_subspaceOverlap<<std::endl;
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
  for (size_t root=0; root<(size_t)m_roots; root++) m_updateShift[root]=-(1+std::numeric_limits<scalar>::epsilon())*m_subspaceEigenvalues[root].real();
}

void Davidson::report()
{
  std::vector<scalar> ev=eigenvalues();
      if (m_verbosity>0) {
        xout << "iteration "<<iterations()<<", error["<<m_worst<<"] = "<<m_error
             << ", eigenvalues: "; for (std::vector<scalar>::const_iterator e=ev.begin(); e!=ev.end(); e++) xout<<" "<<*e;xout<<std::endl;
        }
}


// testing code below here
#include "SimpleParameterVector.h"
static Eigen::MatrixXd testmatrix;

static void _residual(const ParameterVectorSet & psx, ParameterVectorSet & outputs, std::vector<scalar> shift=std::vector<scalar>(), bool append=false) {
  for (size_t k=0; k<psx.size(); k++) {
      Eigen::VectorXd x(testmatrix.rows());
      if (psx[k]->size() != (size_t)testmatrix.rows()) throw std::logic_error("psx wrong size");
      psx[k]->get(&x[0],testmatrix.rows(),0);
      Eigen::VectorXd res = testmatrix * x;
      if (not append) outputs[k]->zero();
      outputs[k]->put(&res[0],testmatrix.rows(),0);
    }
}

static void _preconditoner(const ParameterVectorSet & psg, ParameterVectorSet & psc, std::vector<scalar> shift, bool append=true) {
  size_t n=testmatrix.rows();
  std::vector<scalar> psck(n);
  std::vector<scalar> psgk(n);
  for (size_t k=0; k<psc.size(); k++) {
      psg[k]->get(&psgk[0],n,0);
      if (not append) psc.zero();
      psc[k]->get(&psck[0],n,0);
      for (size_t l=0; l<n; l++)  psck[l] -= psgk[l]/(testmatrix(l,l)+shift[k]);
      psc[k]->put(&psck[0],n,0);
    }
}


typedef SimpleParameterVector ptype;
void Davidson::test(size_t dimension, size_t roots, int verbosity, int problem, bool orthogonalize)
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

  Davidson d(&_residual,&_preconditoner);
  d.m_roots=roots;
  d.m_verbosity=verbosity;
  d.m_maxIterations=dimension;
  d.m_orthogonalize=orthogonalize;
  ParameterVectorSet x;
  ParameterVectorSet g;
  for (size_t root=0; root<(size_t)d.m_roots; root++) {
      ptype* xx=new ptype(dimension);
      xx->zero();
      scalar one=1; xx->put(&one,1,root);
      x.push_back_clone(xx);
      ptype* gg=new ptype(dimension);
      g.push_back_clone(gg);
    }
  xout << "roots="<<roots<<std::endl;
  xout << "initial x="<<x<<std::endl;


  d.solve(g,x);

  std::vector<scalar> ev=d.eigenvalues();
  xout << "Eigenvalues: "; for (std::vector<scalar>::const_iterator e=ev.begin(); e!=ev.end(); e++) xout<<" "<<*e;xout<<std::endl;
  xout << "Reported errors: "; for (std::vector<scalar>::const_iterator e=d.m_errors.begin(); e!=d.m_errors.end(); e++) xout<<" "<<*e;xout<<std::endl;

  _residual(g,x);
  std::vector<scalar> errors;
  for (size_t root=0; root<(size_t)d.m_roots; root++) {
      g[root]->axpy(-ev[root],x[root]);
      errors.push_back(g[root]->dot(g[root]));
    }
  xout << "Square residual norms: "; for (std::vector<scalar>::const_iterator e=errors.begin(); e!=errors.end(); e++) xout<<" "<<*e;xout<<std::endl;
  // be noisy about obvious problems
  if (*std::max_element(errors.begin(),errors.end())>1e-7) throw std::runtime_error("IterativeSolver::Davidson has failed tests");

  for (size_t root=0; root<(size_t)d.m_roots; root++) {
      delete &x[root][0];
      delete &g[root][0];
    }

}
