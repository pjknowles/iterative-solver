#ifndef ITERATIVESOLVER_H
#define ITERATIVESOLVER_H
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <algorithm>
#include <numeric>
#include "LinearAlgebra.h"
#include <Eigen/Dense>
#include <cmath>
#undef isnan
#undef isinf
#ifdef MOLPRO
extern std::ostream &xout;
#else
#define xout std::cout
#endif

namespace LinearAlgebra {
 typedef std::map<std::string, double> optionMap;
 /*!
 * \brief A base class for iterative solvers for linear and non-linear equations, and linear eigensystems.
 *
 * The calling program should set up its own iterative loop, and in each iteration
 * - calculate
 * the action of the matrix on the current expansion vector (linear), or the actual
 * residual (non-linear)
 * - make a call to addVector() which takes the current and previous parameters and proposes
 * an improved estimate, and the best estimate of the residual vector.
 * - calculate a new solution (non-linear) or expansion vector (linear) by implementing
 * appropriate preconditioning on the residual
 * -  make a call to endIteration()
 *
 * Classes that derive from this will, in the simplest case, need to provide just the solveReducedProblem() method that governs how the parameter and residual vectors from successive iterations
 * should be combined to form an optimum solution with minimal residual.
 *
 * The underlying vector spaces are accessed through instances of the vectorSet<scalar> class (or derivatives). These consist of a set of pointers to vector<scalar> objects, which are the vectors themselves; the vectorSet<scalar>
 * object has dimension greater than one in, for example, multi-root diagonalisations where the residual is calculated simultaneously for a number of vectors.
 * Two instances of vectorSet<scalar> have to be passed to iterate() or solve(),
 * and these are used to construct solutions and residuals;
 * this class also creates unlimited additional instances of vectorSet<scalar>,
 * and in memory-sensitive environments, consideration might be given to implementing a derivative of vectorSet<scalar> where the data is resident in external storage.
 * The additional instances are created with hints that they may be stored offline passed to
 * the copy constructor.
 * The use of pointers in vectorSet<scalar> supports polymorphism: the user is free to
 * present instances of classes deriving from vector<scalar>, and thereby implement
 * any special storage arrangements desired.
 */
 template <class scalar=double>
 class IterativeSolver
 {
 public:
  IterativeSolver(
    ) :
    m_PP(0,0),
    m_Pvectors(0),
    m_verbosity(0),
    m_thresh(1e-12),
    m_maxIterations(1000),
    m_minIterations(0),
    m_orthogonalize(false),
    m_linear(false),
    m_hermitian(false),
    m_roots(-1),
    m_date(0),
    m_subspaceMatrixResRes(false),
    m_iterations(0),
    m_singularity_threshold(1e-20)
  {}

  virtual ~IterativeSolver() { }
 public:
  /*!
   * \brief Take, typically, a current solution and residual, and return new solution.
   * In the context of Lanczos-like linear methods, the input will be a current expansion vector and the result of
   * acting on it with the matrix, and the output will be a new expansion vector.
   * For non-linear equations, the input will be the current solution, and the output the current residual.
   * iterate() saves the vectors, calls solveReducedProblem(), calls m_preconditionerFunction(), calls calculateErrors(), and then assesses the error.
   * Derivative classes may often be able to be implemented by changing only solveReducedProblem(), not iterate() or solve().
   * \param parameters On input, the current solution or expansion vector. On exit, the next solution or expansion vector.
   * \param action On input, the residual for solution on entry. On exit, the expected (non-linear) or actual (linear) residual of the interpolated parameters.
   * \param other Optional additional vectors that should be interpolated like the residual.
   * \param options A string of options to be interpreted by solveReducedProblem().
   */
  void addVector( vectorSet<scalar> & parameters, vectorSet<scalar> & action, vectorSet<scalar> & other, const optionMap options=optionMap())
  {
//   xout << "addVector initial parameters: "<<parameters<<std::endl;
   if (m_roots<1) m_roots=parameters.size(); // number of roots defaults to size of parameters
   assert(parameters.size()==action.size());
   m_iterations++;
   for (size_t k=0; k<action.size(); k++) action.m_active[k] = //action.m_active[k] &&
                                                                   parameters.m_active[k];
   m_lastVectorIndex=addVectorSet(parameters,action,other)-1; // derivative classes might eventually store the vectors on top of previous ones, in which case they will need to store the position here for later calculation of iteration step
   solveReducedProblem(parameters,action,other,options);
  }

  void addVector(vectorSet<scalar> & parameters, vectorSet<scalar> & action, const optionMap options=optionMap()) { vectorSet<scalar> other; return addVector(parameters,action,other,options); }

  /*!
   * \brief Specify a P-space vector as a sparse combination of parameters. The container holds a number of segments,
   * each characterised by an offset in the full space, and a vector of coefficients starting at that offset.
   */
  struct Pvector {
   std::vector<size_t> offsets;
   std::vector<std::vector<scalar> > coefficients;
  };

  /*!
   * \brief Add P-space vectors to the expansion set
   * \param Pvectors the vectors to add
   * \param PP Matrix projected onto the existing+new, new P space. It should be provided as a
   * 1-dimensional array, with the existing+new index running fastest.
   */
  void addP(std::vector<Pvector> Pvectors, const scalar* PP ) {
   size_t old=m_PP.rows();
   m_PP.conservativeResize(old+Pvectors.size(),old+Pvectors.size());
   size_t offset=0;
   for (size_t n=0; n<Pvectors.size(); n++) {
    for (int i=0; i<m_PP.rows(); i++)
     m_PP(old+n,i) = m_PP(i,old+n) = PP[offset++];
   }
   for (const auto& Pvector : Pvectors)
    m_Pvectors.push_back(Pvector);
  }

  /*!
     * \brief Take the updated solution vector set, and adjust it if necessary so that it becomes the vector to
     * be used in the next iteration; this is done only in the case of linear solvers where the orthogonalize option is set.
     * Also calculate the degree of convergence, and write progress to std::cout.
     * \param solution The current solution, after interpolation and updating with the preconditioned residual.
     * \param residual The residual after interpolation.
     * \return
     */
  bool endIteration(vectorSet<scalar> & solution, const vectorSet<scalar> & residual)
  {
   calculateErrors(solution,residual);
   adjustUpdate(solution);
   report();
   return m_error < m_thresh;
  }

  /*!
   * \brief Get the solver's suggestion of which degrees of freedom would be best
   * to add to the P-space.
   * \param solution
   * \param residual
   * \param maximumNumber Suggest no more than this number
   * \param threshold Suggest only axes for which the current residual and update
   * indicate an energy improvement in the next iteration of this amount or more.
   * \return
   */
  std::vector<size_t> suggestP(const vectorSet<scalar> & solution,
                               const vectorSet<scalar>& residual,
                               const size_t maximumNumber=1000,
                               const scalar threshold=0)
  {
   return std::vector<size_t>(0);
  }

  /*!
   * \brief Set convergence threshold
   */
  void setThresholds(double thresh) { m_thresh=thresh;}

  unsigned int iterations() { return m_iterations;} //!< How many iterations have occurred

  std::vector<double> eigenvalues() ///< The calculated eigenvalues of m_subspaceMatrix
  {
   std::vector<double> result;
   for (size_t root=0; root<(size_t)m_roots && root < (size_t)m_subspaceEigenvalues.rows(); root++) result.push_back(m_subspaceEigenvalues[root].real());
   return result;
  }

  std::vector<double> errors() {return m_errors;} //!< Error at last iteration

 public:
  Eigen::Matrix<scalar,Eigen::Dynamic,Eigen::Dynamic> m_PP; //!< The PP block of the matrix
  std::vector<Pvector> m_Pvectors;
  int m_verbosity; //!< How much to print.
  double m_thresh; //!< If predicted residual . solution is less than this, converged, irrespective of cthresh and gthresh.
  unsigned int m_maxIterations; //!< Maximum number of iterations in solve()
  unsigned int m_minIterations; //!< Minimum number of iterations in solve()
  bool m_orthogonalize; ///< Whether or not to orthogonalize the result of update() to all previous expansion vectors (appropriate only for linear methods).
  bool m_linear; ///< Whether residuals are linear functions of the corresponding expansion vectors.
  bool m_hermitian; ///< Whether residuals can be assumed to be the action of an underlying self-adjoint operator.
  int m_roots; ///< How many roots to calculate / equations to solve (defaults to size of solution and residual vectors)

 private:
  void adjustUpdate(vectorSet<scalar> & solution)
  {
   //     xout << "m_errors[0] "<<m_errors[0]<<", m_thresh "<<m_thresh<<std::endl;
   for (size_t k=0; k<solution.size(); k++)
    solution.m_active[k] = (m_errors[k] >= m_thresh || m_minIterations>m_iterations);
   //      xout <<  "solution.m_active[0] "<<solution.m_active[0]<<std::endl;
   if (m_orthogonalize) {
    //      xout << "IterativeSolverBase::adjustUpdate solution before orthogonalization: "<<solution<<std::endl;
    for (auto rep=0; rep<2; rep++)
     for (size_t kkk=0; kkk<solution.size(); kkk++) {
      if (solution.m_active[kkk]) {
       for (size_t ll=0; ll<m_solutions.size(); ll++) {
        for (size_t lll=0; lll<m_solutions[ll].size(); lll++) {
         if (m_solutions[ll].m_active[lll]) {
          double s = -(m_solutions[ll][lll]->dot(*solution[kkk])) / (m_solutions[ll][lll]->dot(*m_solutions[ll][lll]));
          solution[kkk]->axpy(s,*m_solutions[ll][lll]);
         }
        }
       }
       for (size_t lll=0; lll<kkk; lll++) {
        if (solution.m_active[lll]) {
         double s = solution[lll]->dot(*solution[kkk]);
         solution[kkk]->axpy(-s,*solution[lll]);
        }
       }
       double s= solution[kkk]->dot(*solution[kkk]);
       if (s <= 0)
        solution.m_active[kkk]=false;
       else
        solution[kkk]->scal(1/std::sqrt(s));
      }
     }
    //      xout << "IterativeSolverBase::adjustUpdate solution after orthogonalization: "<<solution<<std::endl;
   }
  }
 protected:

  virtual void solveReducedProblem(vectorSet<scalar> & residual, vectorSet<scalar> & solution, vectorSet<scalar> & other, const optionMap options=optionMap())=0;
  void solveReducedProblem(vectorSet<scalar> & residual, vectorSet<scalar> & solution, const optionMap options=optionMap()) { vectorSet<scalar> other; solveReducedProblem(residual,solution,other,options); }
  virtual void report()
  {
   if (m_verbosity>0)
    xout << "iteration "<<iterations()<<", error["<<m_worst<<"] = "<<m_error <<std::endl;
  }

 private:
  void calculateSubspaceMatrix(const vectorSet<scalar> &residual, const vectorSet<scalar> &solution)
  {
   //      for (size_t kkk=0; kkk<solution.size(); kkk++)
   //          xout << "solution: "<<solution[kkk]<<std::endl;
   //      for (size_t kkk=0; kkk<residual.size(); kkk++)
   //          xout << "residual: "<<residual[kkk]<<std::endl;
   //      xout << "calculateSubspaceMatrix"<<std::endl;
   //      xout << "residual"<<std::endl<<residual[0]<<std::endl;
   //      xout << "solution"<<std::endl<<solution[0]<<std::endl;
   size_t old_size=m_subspaceMatrix.rows();
   size_t new_size=old_size+std::count(residual.m_active.begin(),residual.m_active.end(),true);
   //      xout << "old_size="<<old_size<<std::endl;
   //      xout << "new_size="<<new_size<<std::endl;
   m_subspaceMatrix.conservativeResize(new_size,new_size);
   m_subspaceOverlap.conservativeResize(new_size,new_size);
   std::vector<vectorSet<scalar>>* bra = m_subspaceMatrixResRes ? &m_residuals : &m_solutions;
   size_t k=old_size;
   for (size_t kkk=0; kkk<residual.size(); kkk++) {
    if (residual.m_active[kkk]) {
     size_t l=0;
     for (size_t ll=0; ll<m_solutions.size(); ll++) {
      for (size_t lll=0; lll<m_solutions[ll].size(); lll++) {
       if (m_solutions[ll].m_active[lll]) {
        //      xout << "bra"<<std::endl<<(*bra)[ll][lll]<<std::endl;
        m_subspaceMatrix(k,l) = m_subspaceMatrix(l,k) = (*bra)[ll][lll]->dot(*residual[kkk]);
        m_subspaceOverlap(k,l) = m_subspaceOverlap(l,k) = m_solutions[ll][lll]->dot(*solution[kkk]);
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

 protected:
  void diagonalizeSubspaceMatrix()
  {
   int kept=m_subspaceMatrix.rows();
   {
    Eigen::EigenSolver<Eigen::Matrix<scalar,Eigen::Dynamic,Eigen::Dynamic> > ss(m_subspaceOverlap);
    Eigen::VectorXcd sse=ss.eigenvalues();
    for (int k=0; k<sse.rows(); k++) {
     if (std::fabs(sse(k).real()) < m_singularity_threshold)
      kept--;
    }
   }
   if (m_verbosity >=0 && kept < m_subspaceMatrix.rows())
    xout <<"IterativeSolver WARNING, subspace singular, pruned from "<<m_subspaceMatrix.rows()<<" to "<<kept<<std::endl;

   Eigen::Matrix<scalar,Eigen::Dynamic,Eigen::Dynamic> H=m_subspaceMatrix.block(0,0,kept,kept);
   Eigen::Matrix<scalar,Eigen::Dynamic,Eigen::Dynamic> S=m_subspaceOverlap.block(0,0,kept,kept);
   Eigen::GeneralizedEigenSolver<Eigen::Matrix<scalar,Eigen::Dynamic,Eigen::Dynamic>> s(H,S);
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
  void calculateErrors(const vectorSet<scalar> & solution, const vectorSet<scalar> &residual)
  {
   if (m_verbosity > 5) {
    xout << "IterativeSolverBase::calculateErrors m_linear" <<m_linear<<std::endl;
    xout << "IterativeSolverBase::calculateErrors solution.m_active"; for (size_t root=0; root<solution.size(); root++) xout <<" "<<solution.m_active[root]; xout <<std::endl;
    xout << "IterativeSolverBase::calculateErrors solution "<<solution<<std::endl;
    xout << "IterativeSolverBase::calculateErrors residual.m_active"; for (size_t root=0; root<solution.size(); root++) xout <<" "<<residual.m_active[root]; xout <<std::endl;
    xout << "IterativeSolverBase::calculateErrors residual "<<residual<<std::endl;
   }
   vectorSet<scalar> step=solution;
   step.axpy(-1,m_solutions[m_lastVectorIndex]);
   if (m_verbosity > 6)
    xout << "IterativeSolverBase::calculateErrors step "<<step<<std::endl;
   m_errors.clear();
   //  xout << "last active "<<m_lastVectorIndex<<" "<<m_residuals[m_lastVectorIndex].m_active[0]<<std::endl;
   for (size_t k=0; k<solution.size(); k++) {
    if (m_linear) // we can use the extrapolated residual if the problem is linear
     m_errors.push_back(residual.m_active[k] ? std::fabs(residual[k]->dot(*step[k])) : 0);
    else
     m_errors.push_back(m_residuals[m_lastVectorIndex].m_active[k] ? std::fabs(m_residuals[m_lastVectorIndex][k]->dot(*step[k])) : 1);
    if (std::isnan(m_errors.back())) throw std::overflow_error("NaN detected in error measure");
   }
   m_error = *max_element(m_errors.begin(),m_errors.end());
   m_worst = max_element(m_errors.begin(),m_errors.end())-m_errors.begin();
   if (m_verbosity > 5) {
    xout << "IterativeSolverBase::calculateErrors m_errors"; for (size_t root=0; root<solution.size(); root++) xout <<" "<<m_errors[root]; xout <<std::endl;
   }
  }

  size_t addVectorSet(const vectorSet<scalar> &solution, const vectorSet<scalar> &residual, const vectorSet<scalar> &other)
  {
   //      if (residual.m_active.front()==0) xout <<"warning: inactive residual"<<std::endl;
   //      if (solution.m_active.front()==0) xout <<"warning: inactive solution"<<std::endl;
   //      for (size_t kkk=0; kkk<solution.size(); kkk++)
   //          xout << "addVectorSet solution: "<<solution[kkk]<<std::endl;
   //      for (size_t kkk=0; kkk<residual.size(); kkk++)
   //          xout << "addVectorSet residual: "<<residual[kkk]<<std::endl;
   m_residuals.emplace_back(vectorSet<scalar>(residual,LINEARALGEBRA_CLONE_ADVISE_DISTRIBUTED|LINEARALGEBRA_CLONE_ADVISE_OFFLINE));
   m_solutions.emplace_back(vectorSet<scalar>(solution,LINEARALGEBRA_CLONE_ADVISE_DISTRIBUTED|LINEARALGEBRA_CLONE_ADVISE_OFFLINE));
   m_others.emplace_back(vectorSet<scalar>(other,LINEARALGEBRA_CLONE_ADVISE_DISTRIBUTED|LINEARALGEBRA_CLONE_ADVISE_OFFLINE));
   calculateSubspaceMatrix(residual,solution);
   return m_residuals.size();
  }

  void deleteVector(size_t index)
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

  std::vector<double> m_errors; //!< Error at last iteration
  double m_error; //!< worst error at last iteration
  size_t m_worst; //!< worst-converged solution, ie m_error = m_errors[m_worst]
  int m_date;
  bool m_subspaceMatrixResRes; // whether m_subspaceMatrix is Residual.Residual (true) or Solution.Residual (false)
  std::vector<vectorSet<scalar>> m_residuals;
  std::vector<vectorSet<scalar>> m_solutions;
  std::vector<vectorSet<scalar>> m_others;
  std::vector<int> m_dateOfBirth;
  size_t m_lastVectorIndex;
  std::vector<scalar> m_updateShift;
  Eigen::Matrix<scalar,Eigen::Dynamic,Eigen::Dynamic> m_subspaceMatrix;
  Eigen::Matrix<scalar,Eigen::Dynamic,Eigen::Dynamic> m_subspaceOverlap;
  Eigen::MatrixXcd m_subspaceEigenvectors; // FIXME templating
  Eigen::VectorXcd m_subspaceEigenvalues; // FIXME templating
 public:
  size_t m_dimension;
 private:
  unsigned int m_iterations;
  double m_singularity_threshold;

 };
}
#include <stdexcept>

namespace LinearAlgebra{

 /** @example LinearEigensystemExample.cpp */
 /*!
 * \brief A class that finds the lowest eigensolutions of a matrix using Davidson's method, i.e. preconditioned Lanczos
 *
 * Example of simplest use: @include LinearEigensystemExample.cpp
 * \tparam scalar Type of matrix elements
 *
 */
 template <class scalar=double>
 class LinearEigensystem : public IterativeSolver<scalar>
 {
 public:
  using IterativeSolver<scalar>::m_verbosity;
  /*!
   * \brief LinearEigensystem
   */
  LinearEigensystem( )
  {
   this->m_linear=true;
   this->m_orthogonalize=true;
  }


 protected:
  virtual void solveReducedProblem(vectorSet<scalar> & solution, vectorSet<scalar> & residual, vectorSet<scalar> & other, const optionMap options=optionMap())
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
          if (m_verbosity>2) xout << "LinearEigensystem::solveReducedProblem kkk="<<kkk<<", ll="<<ll<<", lll="<<lll<<", l="<<l<<std::endl;
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

  virtual void report()
  {
   std::vector<scalar> ev=this->eigenvalues();
   if (m_verbosity>0) {
    xout << "iteration "<<this->iterations()<<", error["<<this->m_worst<<"] = "<<this->m_error
         << ", eigenvalues: "; for (const auto e : ev) xout<<" "<<e;xout<<std::endl;
   }
  }


 };

  /** @example LinearEquationsExample.cpp */
  /*!
  * \brief A class that finds the solutions of linear equation systems using a generalisation of Davidson's method, i.e. preconditioned Lanczos
  *
  * Example of simplest use: @include LinearEquationsExample.cpp
  * \tparam scalar Type of matrix elements
  *
  */
  template <class scalar=double>
  class LinearEquations : public IterativeSolver<scalar>
  {
  public:
   using IterativeSolver<scalar>::m_verbosity;
   /*!
    * \brief Constructor
    * \param rhs right-hand-side vectors. More can be added subsequently using addEquations(), provided iterations have not yet started.
    */
   LinearEquations(const vectorSet<scalar>& rhs = vectorSet<scalar>(0))
   {
    this->m_linear=true;
    this->m_orthogonalize=true;
    addEquations(rhs);
   }

   /*!
    * \brief add one or more equations to the set to be solved, by specifying their right-hand-side vector
    * \param rhs right-hand-side vectors to be added
    */
   void addEquations(const vectorSet<scalar>& rhs) {
    for (const auto& v : rhs ) m_rhs.push_back(v);
   }

  private:
   vectorSet<scalar> m_rhs;


  protected:
   virtual void solveReducedProblem(vectorSet<scalar> & solution, vectorSet<scalar> & residual, vectorSet<scalar> & other, const optionMap options=optionMap())
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
           if (m_verbosity>2) xout << "LinearEquations::solveReducedProblem kkk="<<kkk<<", ll="<<ll<<", lll="<<lll<<", l="<<l<<std::endl;
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

   virtual void report()
   {
    std::vector<scalar> ev=this->eigenvalues();
    if (m_verbosity>0) {
     xout << "iteration "<<this->iterations()<<", error["<<this->m_worst<<"] = "<<this->m_error
          << ", eigenvalues: "; for (const auto e : ev) xout<<" "<<e;xout<<std::endl;
    }
   }


  };

 /** @example DIISexample.cpp */
 /*!
* \brief A class that encapsulates accelerated convergence of non-linear equations
* through the DIIS or related methods.
*
* Example of simplest use: @include DIISexample.cpp
*
*/
template <class scalar=double>
 class DIIS : public IterativeSolver<scalar>
 {
 using IterativeSolver<scalar>::m_residuals;
 using IterativeSolver<scalar>::m_solutions;
 using IterativeSolver<scalar>::m_others;
 public:
 using IterativeSolver<scalar>::m_verbosity;
   enum DIISmode_type {disabled ///< No extrapolation is performed
                       , DIISmode ///< Direct Inversion in the Iterative Subspace
                       , KAINmode ///< Krylov Accelerated Inexact Newton
                      };
   /*!
  * \brief DIIS
  */
 DIIS()
 : m_svdThreshold(1e-10)
 , m_maxDim(6)
{
 this->m_orthogonalize = false;
 setMode(DIISmode);
 Reset();
}

   ~DIIS() { }


   /*!
  * \brief Set options for DIIS.
  * \param mode Whether to perform DIIS, KAIN, or nothing.
  */
   virtual void setMode(enum DIISmode_type mode=DIISmode)
   {
    m_DIISmode = mode;
    this->m_subspaceMatrixResRes = mode!=KAINmode;
//     this->m_preconditionResiduals = mode==KAINmode; // FIXME

    Reset();
    if (m_verbosity>1)
     xout << "m_DIISmode set to "<<m_DIISmode<<std::endl;
   }
   /*!
  * \brief discards previous iteration vectors, but does not clear records
  */
   void Reset()
{
 m_LastResidualNormSq=1e99; // so it can be tested even before extrapolation is done
 this->m_lastVectorIndex=0;
 while (this->m_subspaceMatrix.rows()>0) this->deleteVector(0); ;
 this->m_residuals.clear();
 this->m_solutions.clear();
 this->m_others.clear();
 this->m_Weights.clear();
}
   /*!
  * \brief Introduce a new iteration vector, and perform extrapolation
  * \param residual
  * The vector that
  * will be the one that is analysed to construct the extrapolation.
  * \param solution
  * The current solution that gave rise to residual, and which will be extrapolated to a new predicted solution.
  * \param other (optional)
  * Corresponding other vectors whose sequence will be extrapolated.
  * \param options can contain "weight=xxx" where xxx is the weight to be given to this vector. These options would normally be passed as the corresponding parameter in iterate().
  */
 protected:
   void solveReducedProblem(vectorSet<scalar> & solution, vectorSet<scalar> & residual, vectorSet<scalar> & other, const optionMap options=optionMap())
{
 //	  xout << "Enter DIIS::solveReducedProblem"<<std::endl;
 //	  xout << "residual : "<<residual<<std::endl;
 //	  xout << "solution : "<<solution<<std::endl;
 this->m_updateShift.clear();this->m_updateShift.push_back(-(1+std::numeric_limits<double>::epsilon())*this->m_subspaceMatrix(0,0));
 double weight=options.count("weight") ? (options.find("weight")->second) : 1.0;
 if (this->m_maxDim <= 1 || this->m_DIISmode == disabled) return;

 if (residual.size() > 1) throw std::logic_error("DIIS does not handle multiple solutions");
 this->m_lastVectorIndex=this->m_residuals.size()-1;

 //  if (m_subspaceMatrix.rows() < 9) {
 //      xout << "m_subspaceMatrix on entry to DIIS::solveReducedProblem"<<std::endl<<m_subspaceMatrix<<std::endl;
 //  }
 size_t nDim = this->m_subspaceMatrix.rows();
 this->m_LastResidualNormSq = std::fabs(this->m_subspaceMatrix(nDim-1,nDim-1));
 //  xout << "this->m_LastResidualNormSq "<<this->m_LastResidualNormSq<<std::endl;

 this->m_Weights.push_back(weight);

 Eigen::Array<scalar,Eigen::Dynamic,Eigen::Dynamic> d = this->m_subspaceMatrix.diagonal().array().abs();
 int worst=0, best=0;
 for (size_t i=0; i<nDim; i++) {
     if (d(i) > d(worst)) worst=i;
     if (d(i) < d(best)) best=i;
   }
 double fBaseScale = std::sqrt(d(worst)*d(best));

 while (nDim > this->m_maxDim) { // prune away the worst/oldest vector. Algorithm to be done properly yet
     size_t prune = worst;
     if (true || prune == nDim-1) { // the current vector is the worst, so delete the oldest
         //          xout << "this->m_dateOfBirth: "; for (auto b=this->m_dateOfBirth.begin(); b!=this->m_dateOfBirth.end(); b++) xout <<(*b); xout<<std::endl;
         prune = std::min_element(this->m_dateOfBirth.begin(),this->m_dateOfBirth.end())-this->m_dateOfBirth.begin();
       }
     //      xout << "prune="<<prune<<std::endl;
     //  xout << "nDim="<<nDim<<", m_Weights.size()="<<m_Weights.size()<<std::endl;
     this->deleteVector(prune);
     m_Weights.erase(m_Weights.begin()+prune);
     nDim--;
     //  xout << "nDim="<<nDim<<", m_Weights.size()="<<m_Weights.size()<<std::endl;
   }
 if (nDim != (size_t)this->m_subspaceMatrix.rows()) throw std::logic_error("problem in pruning");
 if (m_Weights.size() != (size_t)this->m_subspaceMatrix.rows()) {
     xout << "nDim="<<nDim<<", m_Weights.size()="<<m_Weights.size()<<std::endl;
     throw std::logic_error("problem after pruning weights");
   }



 //  if (m_subspaceMatrix.rows() < 9) {
 //      xout << "m_subspaceMatrix"<<std::endl<<m_subspaceMatrix<<std::endl;
 //  }

 // build actual DIIS system for the subspace used.
 Eigen::VectorXd
     Rhs(nDim+1),
     Coeffs(nDim+1);
 Eigen::MatrixXd
     B(nDim+1, nDim+1);

 // Factor out common size scales from the residual dots.
 // This is done to increase numerical stability for the case when _all_
 // residuals are very small.
 B.block(0,0,nDim,nDim) = this->m_subspaceMatrix/fBaseScale;
 Rhs.head(nDim) = Eigen::VectorXd::Zero(nDim);

 // make Lagrange/constraint lines.
 for ( size_t i = 0; i < nDim; ++ i )
   B(i, nDim) =  B(nDim, i) = -m_Weights[i];
 B(nDim, nDim) = 0.0;
 Rhs[nDim] = -1;
 //  xout << "B:"<<std::endl<<B<<std::endl;
 //  xout << "Rhs:"<<std::endl<<Rhs<<std::endl;

 // invert the system, determine extrapolation coefficients.
 Eigen::JacobiSVD<Eigen::MatrixXd> svd(B, Eigen::ComputeThinU | Eigen::ComputeThinV);
 svd.setThreshold(m_svdThreshold);
 Coeffs=svd.solve(Rhs).head(nDim);
 m_LastAmplitudeCoeff = Coeffs[nDim-1];
 if (m_verbosity>1) xout << "Combination of iteration vectors: "<<Coeffs.transpose()<<std::endl;
 for (size_t k=0; k<(size_t)Coeffs.rows(); k++)
   if (std::isnan(Coeffs(k))) {
       xout << "B:"<<std::endl<<B<<std::endl;
       xout << "Rhs:"<<std::endl<<Rhs<<std::endl;
       xout << "Combination of iteration vectors: "<<Coeffs.transpose()<<std::endl;
       throw std::overflow_error("NaN detected in DIIS submatrix solution");
     }
 residual.zero();
 solution.zero();
 size_t k=0;
 for (size_t l=0; l<other.size(); l++) other[l]->zero();
 for (size_t kk=0; kk<m_residuals.size(); kk++) {
     if (m_residuals[kk].m_active.front()){
         residual.front()->axpy(Coeffs[k],*m_residuals[kk].front());
         solution.front()->axpy(Coeffs[k],*m_solutions[kk].front());
         for (size_t l=0; l<other.size(); l++)
           other[l]->axpy(Coeffs[k],*m_others[kk][l]);
         k++;
       }
   }
 if (m_verbosity>2) {
     xout << "DIIS.solveReducedProblem() final extrapolated solution: "<<solution.front()<<std::endl;
     xout << "DIIS.solveReducedProblem() final extrapolated residual: "<<residual.front()<<std::endl;
   }
}

   /*!
  * \brief Return the square L2 norm of the extrapolated residual from the last call to solveReducedProblem() or iterate().
  * \return
  */
 public:
   double fLastResidual() const { return m_LastResidualNormSq; }
   /*!
  * \brief Return the coefficient of the last residual vector in the extrapolated residual from the last call to solveReducedProblem() or iterate().
  * \return
  */
   double fLastCoeff() const { return m_LastAmplitudeCoeff; }
   unsigned int nLastDim() const {return m_residuals.size();}
   unsigned int nMaxDim() const { return m_maxDim; }
   double m_svdThreshold; ///< Threshold for singular-value truncation in linear equation solver.
   size_t m_maxDim; ///< Maximum DIIS dimension allowed.
   static void randomTest(size_t sample, size_t n=100, double alpha=0.1, double gamma=0.0, DIISmode_type mode=DIISmode);
 private:
   typedef unsigned int uint;
   enum DIISmode_type m_DIISmode;

   std::vector<double> m_Weights;

   // the following variables are kept for informative/displaying purposes
   double
   // dot(R,R) of last residual vector fed into this state.
   m_LastResidualNormSq,
   // coefficient the actual new vector got in the last DIIS step
   m_LastAmplitudeCoeff;

 };

 /** @example RSPTexample.cpp */
 /*!
* \brief A class that finds the lowest eigensolution of a matrix as a perturbation series
*
* Example of simplest use: @include RSPTexample.cpp
*
*/
template <class scalar>
 class RSPT : public IterativeSolver<scalar>
 {
 using IterativeSolver<scalar>::m_residuals;
 using IterativeSolver<scalar>::m_solutions;
 using IterativeSolver<scalar>::m_others;
 using IterativeSolver<scalar>::m_roots;
 using IterativeSolver<scalar>::m_linear;
 public:
 using IterativeSolver<scalar>::m_verbosity;
   /*!
    * \brief RSPT
    */
 RSPT()
   {
     this->m_linear = true;
     this->m_orthogonalize = false;
     this->m_minIterations=10; // Ensure at least this order
     m_roots=1;
   }
   static void test (size_t n, double alpha);
 protected:
   virtual void solveReducedProblem(vectorSet<scalar> & solution, vectorSet<scalar> & residual, vectorSet<scalar> & other, const optionMap options=optionMap())
{
 size_t n=m_solutions.size();
 // on entry, solution contains |n-1> and residual contains H|n-1>, already stored in m_solutions & m_residuals.
 // on exit, incremental_energies[n] contains E_n = <0 |H-H_0-E_1|n-1> = <0|H|n-1>-(E_0+E_1)<0|n-1>.
 // eventually, Wigner 2n+1 rule should be implemented.
 // on exit, residual contains
 // |d> = -(H_0-E_0)|n>
 //     = (H-H_0-E_1)|n-1> - sum_{k=0}^{n-2} E_{n-k} |k>.
 //     = (H-E_0-E_1)|n-1> - (H_0-E_0)|n-1> - sum_{k=0}^{n-2} E_{n-k} |k>.
 // on exit, solution contains 0.

 {
 Eigen::MatrixXcd oldeigenvalues=this->m_subspaceEigenvalues;
 this->diagonalizeSubspaceMatrix();
 if (std::isnan(this->m_subspaceEigenvalues(0).real())) this->m_subspaceEigenvalues=oldeigenvalues;
 }
 if (m_verbosity>1) xout << "Subspace eigenvalues"<<std::endl<<this->m_subspaceEigenvalues<<std::endl;
 if (m_verbosity>2) xout << "Subspace eigenvectors"<<std::endl<<this->m_subspaceEigenvectors<<std::endl;

//      xout << "|"<<n-1<<">: "<<solution;
//      xout << "H|"<<n-1<<">: "<<residual;
 if (n == 1) {
     // the preconditioner function must take special action when shift=0
     // and return H0.residual not -(H0-shift)^{-1}.residual
     std::vector<double> shift; shift.push_back(0);
     throw std::logic_error("RSPT coding not complete");
//      m_preconditionerFunction(solution,residual,shift,false);
//      xout << "solution="<<solution;
//      xout << "residual="<<residual;
     m_E0 = solution.front()->dot(*residual.front());
 }
     solution.zero();
     residual.zero();
     residual.axpy(1,m_residuals.back());
 if (n == 1) {
     m_incremental_energies.resize(2);
     m_incremental_energies[0]=m_E0;
     m_incremental_energies[1]=this->m_subspaceMatrix(0,0)-m_incremental_energies[0];
     residual.axpy(-this->m_subspaceMatrix(0,0),m_solutions.back());
     this->m_updateShift.clear();this->m_updateShift.push_back(-(1+std::numeric_limits<double>::epsilon())*m_E0);
     if (m_verbosity>=0)
         xout << "E_0="<<m_E0<<std::endl << "E_1="<<m_incremental_energies[1]<<std::endl;
//      xout << "d(1)after incrementing solution"<<residual<<std::endl;
 } else {
     m_incremental_energies.push_back(this->m_subspaceMatrix(n-1,0)-this->m_subspaceMatrix(0,0)*this->m_subspaceOverlap(n-1,0));
     if (m_verbosity>=0)
         xout << "E_"<<n<<"="<<m_incremental_energies.back()<<", Eigenvalue="<<this->m_subspaceEigenvalues(0)<<std::endl;
//      xout << "d(n)=g(n-1)"<<residual<<std::endl;
     residual.axpy(1,m_lastH0mE0psi);
//      xout << "d(n)=g(n-1)+d(n-1)"<<residual<<std::endl;
//      xout << "H(0,0) "<<this->m_subspaceMatrix(0,0)<<std::endl;
//      xout << "c(n-1) "<<m_solutions.back()<<std::endl;
     residual.axpy(-this->m_subspaceMatrix(0,0),m_solutions.back());
//      xout << "d(n)=g(n-1)+d(n-1)-(E0+E1)c(n-1)"<<residual<<std::endl;
 // this is structured for multistate, but not thought about yet
 for (size_t kkk=0; kkk<residual.size(); kkk++) {
     size_t l=0;
     for (size_t ll=0; ll<n-1; ll++) {
         for (size_t lll=0; lll<m_solutions[ll].size(); lll++) {
                 residual[kkk]->axpy(-m_incremental_energies[n-ll],*m_solutions[ll][lll]);
//      xout << "d(n)after incrementing solution"<<residual<<std::endl;
                 l++;
           }
       }
   }

 }
 m_lastH0mE0psi = residual; // we will need this in the next iteration // FIXME does this leak memory?
//  xout << "-(H0-E0)|"<<n<<">: "<<residual<<std::endl;
}

 public:
   int m_order; ///< Up to what order of perturbation theory should the energy be obtained.
   std::vector<double> incremental_energies(size_t state=0) {return m_incremental_energies;} ///< The incremental energies order by order.
///> The total energy to a given order.
   double energy(size_t order, ///< the desired maximum order of perturbation theory
                 size_t state=0 ///< the desired state
                              )
{
   double result=0.0;
   for (size_t k=0; k<=order; k++)
       result += m_incremental_energies[k];
   return result;
}

 private:
   vectorSet<scalar> m_lastH0mE0psi;
   std::vector<double> m_incremental_energies;
   double m_E0;
 };

}



  // C interface
 extern "C" void IterativeSolverLinearEigensystemInitialize(size_t nQ, size_t nroot, double thresh, int maxIterations, int verbosity);
 extern "C" void IterativeSolverLinearEigensystemAddVector(double* c, double* g, double* eigenvalue);
 extern "C" int IterativeSolverLinearEigensystemEndIteration(double* c, double* g, double* error);
 extern "C" void IterativeSolverAddP(size_t* indices, double* coefficients, double* pp);



#endif // ITERATIVESOLVER_H
