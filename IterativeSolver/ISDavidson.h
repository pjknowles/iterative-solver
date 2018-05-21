#ifndef DAVIDSON_H
#define DAVIDSON_H
#include "IterativeSolver.h"
#include <stdexcept>
#include "PagedVector.h"

namespace LinearAlgebra{

 /** @example DavidsonExample.cpp */
 /*!
 * \brief A class that finds the lowest eigensolutions of a matrix using Davidson's method, i.e. preconditioned Lanczos
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
          if (m_verbosity>2) xout << "Davidson::extrapolate kkk="<<kkk<<", ll="<<ll<<", lll="<<lll<<", l="<<l<<std::endl;
          if (m_verbosity>2) xout << "Eigenvectors:\n"<<this->m_subspaceEigenvectors(l,kkk).real()<<std::endl;
       solution[kkk]->axpy(this->m_subspaceEigenvectors(l,kkk).real(),*this->m_solutions[ll][lll]);
       residual[kkk]->axpy(this->m_subspaceEigenvectors(l,kkk).real(),*this->m_residuals[ll][lll]);
       l++;
      }
     }
    }
    residual[kkk]->axpy(-this->m_subspaceEigenvalues(kkk).real(),*solution[kkk]);
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
  // C interface
  using v = PagedVector<double>;

  static std::unique_ptr<Davidson<double> > instance;
  extern "C" void IterativeSolverDavidsonInitialize(size_t nQ, size_t nP, size_t nroot, double* PP) {
   instance.reset(new Davidson<double>(Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> >(PP,nP,nP)));
   instance->m_dimension=nQ;
   instance->m_roots=nroot;
  }
  extern "C" void IterativeSolverDavidsonInterpolate(double* c, double* g, double* eigenvalue) {
   vectorSet<double> cc,gg;
   for (int root=0; root < instance->m_roots; root++) {
    cc.push_back(std::shared_ptr<v>(new v(instance->m_dimension)));
    cc.back()->put(&c[root*instance->m_dimension],instance->m_dimension,0);
    gg.push_back(std::shared_ptr<v>(new v(instance->m_dimension)));
    gg.back()->put(&g[root*instance->m_dimension],instance->m_dimension,0);
   }
   instance->addVector(cc,gg);
   for (int root=0; root < instance->m_roots; root++) {
    cc[root]->get(&c[root*instance->m_dimension],instance->m_dimension,0);
    gg[root]->get(&g[root*instance->m_dimension],instance->m_dimension,0);
    eigenvalue[root] = instance->eigenvalues()[root];
   }
  }

  extern "C" int IterativeSolverDavidsonFinalize(double* c, double* g, double* error) {
   vectorSet<double> cc,gg;
   for (int root=0; root < instance->m_roots; root++) {
    cc.push_back(std::shared_ptr<v>(new v(instance->m_dimension)));
    cc.back()->put(&c[root*instance->m_dimension],instance->m_dimension,0);
    gg.push_back(std::shared_ptr<v>(new v(instance->m_dimension)));
    gg.back()->put(&g[root*instance->m_dimension],instance->m_dimension,0);
   }
   bool result = instance->finalize(cc,gg);
   for (int root=0; root < instance->m_roots; root++) {
    cc[root]->get(&c[root*instance->m_dimension],instance->m_dimension,0);
    error[root] = instance->errors()[root];
   }
   return result;
  }

}

#endif // DAVIDSON_H
