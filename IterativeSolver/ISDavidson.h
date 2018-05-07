#ifndef DAVIDSON_H
#define DAVIDSON_H
#include "IterativeSolver.h"
#include <stdexcept>

namespace LinearAlgebra{

 /** @example DavidsonExample.cpp */
 /*!
 * \brief A class that finds the lowest eigensolutions of a matrix using Davidson's method
 *
 * Example of simplest use: @include DavidsonExample.cpp
 * \tparam scalar Type of matrix elements
 *
 */
 template <class scalar=double>
 class Davidson : public IterativeSolverBase<scalar>
 {
 public:
  using IterativeSolverBase<scalar>::m_verbosity;
  /*!
   * \brief Davidson
   * \param PP The PP block of the matrix
   */
  Davidson( const Eigen::Matrix<scalar,Eigen::Dynamic,Eigen::Dynamic>& PP=Eigen::Matrix<scalar,Eigen::Dynamic,Eigen::Dynamic>(0,0) )
   : IterativeSolverBase<scalar>(PP)
  {
   this->m_linear=true;
   this->m_orthogonalize=true;
  }


 protected:
  virtual void extrapolate(vectorSet<scalar> & solution, vectorSet<scalar> & residual, vectorSet<scalar> & other, const optionMap options=optionMap())
  {
   if (m_verbosity>2) xout << "Subspace matrix"<<std::endl<<this->m_subspaceMatrix<<std::endl;
   if (m_verbosity>2) xout << "Subspace overlap"<<std::endl<<this->m_subspaceOverlap<<std::endl;
   this->diagonalizeSubspaceMatrix();

   if (m_verbosity>1) xout << "Subspace eigenvalues"<<std::endl<<this->m_subspaceEigenvalues<<std::endl;
   if (m_verbosity>2) xout << "Subspace eigenvectors"<<std::endl<<this->m_subspaceEigenvectors<<std::endl;
   residual.zero();
   solution.zero();
   for (size_t kkk=0; kkk<residual.size(); kkk++) {
    size_t l=0;
    for (size_t ll=0; ll<this->m_solutions.size(); ll++) {
     for (size_t lll=0; lll<this->m_solutions[ll].size(); lll++) {
      if (this->m_solutions[ll].m_active[lll]) {
       solution[kkk]->axpy(this->m_subspaceEigenvectors(l,kkk).real(),this->m_solutions[ll][lll]);
       residual[kkk]->axpy(this->m_subspaceEigenvectors(l,kkk).real(),this->m_residuals[ll][lll]);
       l++;
      }
     }
    }
    residual[kkk]->axpy(-this->m_subspaceEigenvalues(kkk).real(),solution[kkk]);
   }

   this->m_updateShift.resize(this->m_roots);
   for (size_t root=0; root<(size_t)this->m_roots; root++) this->m_updateShift[root]=-(1+std::numeric_limits<scalar>::epsilon())*this->m_subspaceEigenvalues[root].real();
  }

  virtual void extrapolate(vectorSet<scalar> & residual, vectorSet<scalar> & solution, const optionMap options=optionMap()) { vectorSet<scalar> other; extrapolate(solution,residual,other,options); }

  virtual void report()
  {
   std::vector<scalar> ev=this->eigenvalues();
   if (m_verbosity>0) {
    xout << "iteration "<<this->iterations()<<", error["<<this->m_worst<<"] = "<<this->m_error
         << ", eigenvalues: "; for (const auto e : ev) xout<<" "<<e;xout<<std::endl;
   }
  }


 };
  // testing code below here
  /*!
   * \brief Test
   * \param dimension The dimension of the test matrix
   * \param roots How many eigensolutions to find
   * \param verbosity How much to report
   * \param problem Selects which test matrix to use
   * \param orthogonalize Whether to orthogonalize expansion vectors
   * \tparam ptype Concrete class template that implements LinearAlgebra::vectorSet
   * \tparam scalar Type of matrix elements
   */
 template <class ptype, class scalar=double>
  static void DavidsonTest(size_t dimension, size_t roots=1, int verbosity=0, int problem=0, bool orthogonalize=true)
  {

   static Eigen::Matrix<scalar,Eigen::Dynamic,Eigen::Dynamic> testmatrix;

   static struct {
    void operator()(const vectorSet<scalar> & psx, vectorSet<scalar> & outputs) const {
     for (size_t k=0; k<psx.size(); k++) {
      Eigen::VectorXd x(testmatrix.rows());
      if (psx[k]->size() != (size_t)testmatrix.rows()) throw std::logic_error("psx wrong size");
      psx[k]->get(&x[0],testmatrix.rows(),0);
      Eigen::VectorXd res = testmatrix * x;
      outputs[k]->put(&res[0],testmatrix.rows(),0);
     }
    }
   } action;

   static struct {
    void operator()(vectorSet<scalar> & psc, const vectorSet<scalar> & psg, std::vector<scalar> shift, bool append=true) const {
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
   } update;

   xout << "Test IterativeSolver::Davidson dimension="<<dimension<<", roots="<<roots<<", problem="<<problem<<", orthogonalize="<<orthogonalize<<std::endl;
   testmatrix.resize(dimension,dimension);
   for (size_t k=0; k<dimension; k++)
    for (size_t l=0; l<dimension; l++)
     if (problem==0)
      testmatrix(l,k)=-1;
     else if (problem==1)
      testmatrix(l,k)=l+k+2;
     else if (problem==2)
      testmatrix(l,k)=( k==l ? k+1 : 1);
     else if (problem==3)
      testmatrix(l,k)=( k==l ? 1 : 1);
     else
      throw std::logic_error("invalid problem in Davidson::test");
   if (problem==3) testmatrix(0,1)=testmatrix(1,0)=1;

   Davidson<scalar> d;
   d.m_roots=roots;
   d.m_verbosity=verbosity;
   d.m_maxIterations=dimension;
   d.m_orthogonalize=orthogonalize;
   vectorSet<scalar> x;
   vectorSet<scalar> g;
   for (size_t root=0; root<(size_t)d.m_roots; root++) {
    auto xx=std::make_shared<ptype>(dimension);
    xx->zero();
    scalar one=1; xx->put(&one,1,root);
    x.push_back_clone(xx);
    auto gg=std::make_shared<ptype>(dimension);
    g.push_back_clone(gg);
   }

   for (size_t iteration=0; iteration<dimension; iteration++) {
    action(x,g);
    d.interpolate(x,g);
    std::vector<scalar> shift;
    for (size_t root=0; root<(size_t)d.m_roots; root++) shift.push_back(-d.eigenvalues()[root]+1e-14);
    update(x,g,shift);
    if (d.finalize(x,g)) break;
   }
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<scalar,Eigen::Dynamic,Eigen::Dynamic>> es(testmatrix);
    xout << "true eigenvalues: "<<es.eigenvalues().head(d.m_roots).transpose()<<std::endl;
//    xout << "true eigenvectors:\n"<<es.eigenvectors().leftCols(d.m_roots).transpose()<<std::endl;


   std::vector<scalar> ev=d.eigenvalues();
   xout << "Eigenvalues: "; size_t root=0; for (const auto& e : ev) xout<<" "<<e<<"(error="<<e-es.eigenvalues()(root++)<<")";xout<<std::endl;
   xout << "Reported errors: "; for (const auto& e: d.errors()) xout<<" "<<e;xout<<std::endl;

   action(x,g);
   std::vector<scalar> errors;
   for (size_t root=0; root<(size_t)d.m_roots; root++) {
    g[root]->axpy(-ev[root],x[root]);
    errors.push_back(g[root]->dot(g[root]));
   }
//   xout << "Square residual norms: "; for (typename std::vector<scalar>::const_iterator e=errors.begin(); e!=errors.end(); e++) xout<<" "<<*e;xout<<std::endl;
   xout << "Square residual norms: "; for (const auto& e: errors) xout<<" "<<e;xout<<std::endl;
   // be noisy about obvious problems
   if (*std::max_element(errors.begin(),errors.end())>1e-7) throw std::runtime_error("IterativeSolver::Davidson has failed tests");

  }

}

#endif // DAVIDSON_H
