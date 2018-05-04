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

  /*!
 * \brief A base class for iterative solvers such as DIIS, KAIN, Davidson. The class provides support for preconditioned update, via a provided functor.
 *
 * The user needs to provide the two routines residualFunction() and preconditionerFunction() through the class constructor. These define the problem being solved: the first should calculate the residual
 * or action vector from a solution vector,
 * and the second should apply the negative of a preconditioner to a provided residual vector, optionally adding it to existing contents of a result vector.  The user also needs to provide an initial guess in the call to solve() or iterate().
 *
 * Two drivers are provided: the calling program can set up its own iterative loop, and in each loop call residualFunction() and iterate(); this gives the flexibility to pass additional parameters
 * to residualFunction(). The simpler mode of use is a single call to solve(), which manages the iterations itself.
 *
 * Classes that derive from this will, in the simplest case, need to provide just the extrapolate() method that governs how the solution and residual vectors from successive iterations
 * should be combined to form an optimum solution with minimal residual.  In more complicated cases - for example, in Davidson's method, where the preconditioner depends on the current energy -
 * it will be necessary to reimplement also the iterate() method.
 *
 * The underlying vector spaces are accessed through instances of the ParameterVectorSet class (or derivatives). These consist of a set of ParameterVector objects, which are the vectors themselves; the ParameterVectorSet
 * object has dimension greater than one in, for example, multi-root diagonalisations where the residual is calculated simultaneously for a number of vectors.
 * Two instances of ParameterVectorSet have to be passed to iterate() or solve(), and these are used to construct solutions and residuals; this class also creates unlimited additional instances of ParameterVectorSet,
 * and in memory-sensitive environments, consideration might be given to implementing a derivative of ParameterVectorSet where the data is resident in external storage.
 */
template <class scalar=double>
  class IterativeSolverBase
  {
  public:
  typedef vector<scalar> ParameterVector;
  typedef vectorSet<scalar> ParameterVectorSet;
    //    typedef void (*ParameterSetTransformation)(const ParameterVectorSet & inputs, ParameterVectorSet & outputs, std::vector<scalar> shift, bool append);
    struct ParameterSetTransformation
    {
      virtual void operator()(const ParameterVectorSet & inputs, ParameterVectorSet & outputs, std::vector<scalar> shift, bool append) const =0;
    };
    /*!
 * \brief Place-holding template for residual calculation. It just returns the input as output.
 */
    static struct : ParameterSetTransformation {
     /*!
       * \brief operator()
 * \param inputs The parameters.
 * \param outputs On output, contains the corresponding residual vectors.
 * \param shift
 * \param append Whether to add the result to the original content of outputs
       */
      void operator()( const ParameterVectorSet & inputs, ParameterVectorSet & outputs, std::vector<scalar> shift=std::vector<scalar>(), bool append=false) const override { if (append) outputs=inputs; }
    } noOp ;
    /*!
 * \brief Place-holding template for update calculation. It just returns the input as output.
 */
    static struct : ParameterSetTransformation {
     /*!
       * \brief operator()
 * \param inputs The parameters.
 * \param append Whether to add the result to the original content of outputs
 * \param outputs On output, contains the corresponding residual vectors.
 * \param shift
       */
      void operator()(const ParameterVectorSet & inputs, ParameterVectorSet & outputs, std::vector<scalar> shift=std::vector<scalar>(), bool append=true) const override {
        if (not append) outputs.zero();
        for (size_t k=0; k<inputs.size(); k++)
          outputs[k]->axpy(-1,inputs[k]);
      }
    } steepestDescent ;

    /*!
   * \brief IterativeSolverBase
   * \param residualFunction A functor that evaluates the residual vectors. Used by method solve(); does not have to be provided if iterations are constructed explicitly in the calling program.
   * \param preconditionerFunction A functor that applies a preconditioner to a residual to give an update. Used by methods iterate() and solve().
   */
    IterativeSolverBase(const ParameterSetTransformation& residualFunction, const ParameterSetTransformation& preconditionerFunction=steepestDescent)
      :  m_preconditionerFunction(preconditionerFunction),
        m_residualFunction(residualFunction),
        m_verbosity(0),
        m_thresh(1e-12),
        m_maxIterations(1000),
        m_minIterations(0),
        m_orthogonalize(false),
        m_linear(false),
        m_hermitian(false),
        m_preconditionResiduals(false),
        m_roots(-1),
        m_date(0),
        m_subspaceMatrixResRes(false),
        m_iterations(0),
        m_singularity_threshold(1e-20)
    {}

    virtual ~IterativeSolverBase() { }
  public:
    /*!
   * \brief The functor that will take the current solution and residual, and produce the predicted solution.
   */
    const ParameterSetTransformation& m_preconditionerFunction;
    /*!
   * \brief The functor that will take a current solution and calculate the residual.
   */
    const ParameterSetTransformation& m_residualFunction;
  public:
    /*!
   * \brief Take, typically, a current solution and residual, and return new solution.
   * In the context of Lanczos-like methods, the input will be a current expansion vector and the result of
   * acting on it with the matrix, and the output will be a new expansion vector.
   * iterate() saves the vectors, calls extrapolate(), calls m_preconditionerFunction(), calls calculateErrors(), and then assesses the error.
   * Derivative classes may often be able to be implemented by changing only extrapolate(), not iterate() or solve().
   * \param residual On input, the residual for solution on entry. On exit, the extrapolated residual.
   * \param solution On input, the current solution or expansion vector. On exit, the next solution or expansion vector.
   * \param other Optional additional vectors that should be extrapolated.
   * \param options A string of options to be interpreted by extrapolate().
   */
    typedef std::map<std::string, double> optionMap;
    virtual bool iterate(ParameterVectorSet & residual, ParameterVectorSet & solution, ParameterVectorSet & other, const optionMap options=optionMap())
    {
      if (m_roots<1) m_roots=solution.size(); // number of roots defaults to size of solution
      assert(solution.size()==residual.size());
      m_iterations++;
      for (size_t k=0; k<residual.size(); k++) residual.m_active[k] = residual.m_active[k] && solution.m_active[k];
      if (m_preconditionResiduals) m_preconditionerFunction(residual,residual,m_updateShift,false);
      m_lastVectorIndex=addVectorSet(residual,solution,other)-1; // derivative classes might eventually store the vectors on top of previous ones, in which case they will need to store the position here for later calculation of iteration step
      extrapolate(residual,solution,other,options);
      if (m_preconditionResiduals)
        solution.axpy(1,residual);
      else
        m_preconditionerFunction(residual,solution,m_updateShift,true);
      calculateErrors(solution,residual);
      adjustUpdate(solution);
      return m_error < m_thresh;
    }

    virtual bool iterate(ParameterVectorSet & residual, ParameterVectorSet & solution, const optionMap options=optionMap()) { ParameterVectorSet other; return iterate(residual,solution,other,options); }
    /*!
   * \brief Solve iteratively by repeated calls to residualFunction() and iterate().
   * \param residual Ignored on input; on exit, contains the final residual; used as working space.
   * \param solution On input, contains an initial guess; on exit, contains the final solution.
   * \param options A string of options to be interpreted by extrapolate().
   * \return Whether or not convergence has been reached.
   */
    virtual bool solve(ParameterVectorSet & residual, ParameterVectorSet & solution, const optionMap options=optionMap());
    /*!
   * \brief Set convergence threshold
   */


    void setThresholds(double thresh) { m_thresh=thresh;}

    unsigned int iterations() { return m_iterations;} //!< How many iterations have occurred
    std::vector<double> eigenvalues(); ///< The calculated eigenvalues of m_subspaceMatrix
    std::vector<double> errors() {return m_errors;} //!< Error at last iteration

  public:
    int m_verbosity; //!< How much to print.
    double m_thresh; //!< If predicted residual . solution is less than this, converged, irrespective of cthresh and gthresh.
    unsigned int m_maxIterations; //!< Maximum number of iterations in solve()
    unsigned int m_minIterations; //!< Minimum number of iterations in solve()
    bool m_orthogonalize; ///< Whether or not to orthogonalize the result of update() to all previous expansion vectors (appropriate only for linear methods).
    bool m_linear; ///< Whether residuals are linear functions of the corresponding expansion vectors.
    bool m_hermitian; ///< Whether residuals can be assumed to be the action of an underlying self-adjoint operator.
    bool m_preconditionResiduals; ///< Whether the subspace algorithm should work with preconditioned or raw residual vectors
    int m_roots; ///< How many roots to calculate / equations to solve (defaults to size of solution and residual vectors)

  protected:
    virtual void adjustUpdate(ParameterVectorSet & solution);
    virtual void extrapolate(ParameterVectorSet & residual, ParameterVectorSet & solution, ParameterVectorSet & other, const optionMap options=optionMap())=0;
    virtual void extrapolate(ParameterVectorSet & residual, ParameterVectorSet & solution, const optionMap options=optionMap()) { ParameterVectorSet other; extrapolate(residual,solution,other,options); }
    virtual void report();
    void calculateSubspaceMatrix(const ParameterVectorSet &residual, const ParameterVectorSet &solution);
    void diagonalizeSubspaceMatrix();
    void calculateErrors(const ParameterVectorSet & solution, const ParameterVectorSet &residual);
    size_t addVectorSet(const ParameterVectorSet &residual, const ParameterVectorSet &solution, const ParameterVectorSet &other);
    void deleteVector(size_t index);
    std::vector<double> m_errors; //!< Error at last iteration
    double m_error; //!< worst error at last iteration
    size_t m_worst; //!< worst-converged solution, ie m_error = m_errors[m_worst]
    int m_date;
    bool m_subspaceMatrixResRes; // whether m_subspaceMatrix is Residual.Residual (true) or Solution.Residual (false)
    std::vector<ParameterVectorSet> m_residuals;
    std::vector<ParameterVectorSet> m_solutions;
    std::vector<ParameterVectorSet> m_others;
    std::vector<int> m_dateOfBirth;
    size_t m_lastVectorIndex;
    std::vector<scalar> m_updateShift;
    Eigen::MatrixXd m_subspaceMatrix;
    Eigen::MatrixXd m_subspaceOverlap;
    Eigen::MatrixXcd m_subspaceEigenvectors;
    Eigen::VectorXcd m_subspaceEigenvalues;
  private:
    unsigned int m_iterations;
    double m_singularity_threshold;







    bool IterativeSolverBase::solve(ParameterVectorSet & residual, ParameterVectorSet & solution, const optionMap options)
    {
      bool converged=false;
      for (unsigned int iteration=1; iteration <= m_maxIterations && (not converged || iteration <= m_minIterations); iteration++) {
          m_residualFunction(solution,residual,std::vector<double>(),false);
          converged = iterate(residual,solution);
          report();
        }
      return converged;
    }

    void IterativeSolverBase::report()
    {
          if (m_verbosity>0)
            xout << "iteration "<<iterations()<<", error["<<m_worst<<"] = "<<m_error <<std::endl;
    }

    void IterativeSolverBase::adjustUpdate(ParameterVectorSet &solution)
    {
      for (size_t k=0; k<solution.size(); k++)
        solution.m_active[k] = (m_errors[k] >= m_thresh || m_minIterations>m_iterations);
      if (m_orthogonalize) {
          //      xout << "IterativeSolverBase::adjustUpdate solution before orthogonalization: "<<solution<<std::endl;
          for (auto rep=0; rep<2; rep++)
          for (size_t kkk=0; kkk<solution.size(); kkk++) {
              if (solution.m_active[kkk]) {
                  for (size_t ll=0; ll<m_solutions.size(); ll++) {
                      for (size_t lll=0; lll<m_solutions[ll].size(); lll++) {
                          if (m_solutions[ll].m_active[lll]) {
                              double s = -(m_solutions[ll][lll]->dot(solution[kkk])) / (m_solutions[ll][lll]->dot(m_solutions[ll][lll]));
                              solution[kkk]->axpy(s,m_solutions[ll][lll]);
                            }
                        }
                    }
                  for (size_t lll=0; lll<kkk; lll++) {
                      if (solution.m_active[lll]) {
                          double s = solution[lll]->dot(solution[kkk]);
                          solution[kkk]->axpy(-s,solution[lll]);
                        }
                    }
                  double s= solution[kkk]->dot(solution[kkk]);
                  if (s <= 0)
                    solution.m_active[kkk]=false;
                  else
                    solution[kkk]->scal(1/std::sqrt(s));
                }
            }
          //      xout << "IterativeSolverBase::adjustUpdate solution after orthogonalization: "<<solution<<std::endl;
        }
    }

    size_t IterativeSolverBase::addVectorSet(const ParameterVectorSet &residual, const ParameterVectorSet &solution, const ParameterVectorSet &other)
    {
    //  if (residual.m_active.front()==0) xout <<"warning: inactive residual"<<std::endl;
    //  if (solution.m_active.front()==0) xout <<"warning: inactive solution"<<std::endl;
        m_residuals.push_back(residual);
        m_solutions.push_back(solution);
    //    m_residuals.emplace_back(new ParameterVector(residual,LINEARALGEBRA_CLONE_ADVISE_DISTRIBUTED|LINEARALGEBRA_CLONE_ADVISE_OFFLINE));
    //    m_solutions.emplace_back(new ParameterVector(solution,LINEARALGEBRA_CLONE_ADVISE_DISTRIBUTED|LINEARALGEBRA_CLONE_ADVISE_OFFLINE));
        m_others.push_back(other);
        calculateSubspaceMatrix(residual,solution);
        return m_residuals.size();
    }

    void IterativeSolverBase::deleteVector(size_t index)
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

    void IterativeSolverBase::calculateSubspaceMatrix(const ParameterVectorSet &residual, const ParameterVectorSet &solution)
    {
    //  xout << "calculateSubspaceMatrix"<<std::endl;
    //  xout << "residual"<<std::endl<<residual[0]<<std::endl;
    //  xout << "solution"<<std::endl<<solution[0]<<std::endl;
      size_t old_size=m_subspaceMatrix.rows();
      size_t new_size=old_size+std::count(residual.m_active.begin(),residual.m_active.end(),true);
      m_subspaceMatrix.conservativeResize(new_size,new_size);
      m_subspaceOverlap.conservativeResize(new_size,new_size);
      std::vector<ParameterVectorSet>* bra = m_subspaceMatrixResRes ? &m_residuals : &m_solutions;
      size_t k=old_size;
      for (size_t kkk=0; kkk<residual.size(); kkk++) {
          if (residual.m_active[kkk]) {
              size_t l=0;
              for (size_t ll=0; ll<m_solutions.size(); ll++) {
                  for (size_t lll=0; lll<m_solutions[ll].size(); lll++) {
                      if (m_solutions[ll].m_active[lll]) {
    //  xout << "bra"<<std::endl<<(*bra)[ll][lll]<<std::endl;
                          m_subspaceMatrix(k,l) = m_subspaceMatrix(l,k) = (*bra)[ll][lll]->dot(residual[kkk]);
                          m_subspaceOverlap(k,l) = m_subspaceOverlap(l,k) = m_solutions[ll][lll]->dot(solution[kkk]);
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

    void IterativeSolverBase::diagonalizeSubspaceMatrix()
    {
      int kept=m_subspaceMatrix.rows();
      {
          Eigen::EigenSolver<Eigen::MatrixXd> ss(m_subspaceOverlap);
          Eigen::VectorXcd sse=ss.eigenvalues();
          for (int k=0; k<sse.rows(); k++) {
              if (std::fabs(sse(k).real()) < m_singularity_threshold)
                  kept--;
          }
      }
      if (m_verbosity >=0 && kept < m_subspaceMatrix.rows())
          xout <<"IterativeSolver WARNING, subspace singular, pruned from "<<m_subspaceMatrix.rows()<<" to "<<kept<<std::endl;

      Eigen::MatrixXd H=m_subspaceMatrix.block(0,0,kept,kept);
      Eigen::MatrixXd S=m_subspaceOverlap.block(0,0,kept,kept);
      Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> s(H,S);
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


    void IterativeSolverBase::calculateErrors(const ParameterVectorSet &solution, const ParameterVectorSet &residual)
    {
      if (m_verbosity > 5) {
        xout << "IterativeSolverBase::calculateErrors m_linear" <<m_linear<<std::endl;
        xout << "IterativeSolverBase::calculateErrors solution.m_active"; for (size_t root=0; root<solution.size(); root++) xout <<" "<<solution.m_active[root]; xout <<std::endl;
        xout << "IterativeSolverBase::calculateErrors solution "<<solution<<std::endl;
        xout << "IterativeSolverBase::calculateErrors residual.m_active"; for (size_t root=0; root<solution.size(); root++) xout <<" "<<residual.m_active[root]; xout <<std::endl;
        xout << "IterativeSolverBase::calculateErrors residual "<<residual<<std::endl;
        }
      ParameterVectorSet step=solution;
      step.axpy(-1,m_solutions[m_lastVectorIndex]);
      if (m_verbosity > 6)
        xout << "IterativeSolverBase::calculateErrors step "<<step<<std::endl;
      m_errors.clear();
    //  xout << "last active "<<m_lastVectorIndex<<" "<<m_residuals[m_lastVectorIndex].m_active[0]<<std::endl;
      for (size_t k=0; k<solution.size(); k++) {
        if (m_linear) // we can use the extrapolated residual if the problem is linear
          m_errors.push_back(residual.m_active[k] ? std::fabs(residual[k]->dot(step[k])) : 0);
        else
          m_errors.push_back(m_residuals[m_lastVectorIndex].m_active[k] ? std::fabs(m_residuals[m_lastVectorIndex][k]->dot(step[k])) : 1);
        if (std::isnan(m_errors.back())) throw std::overflow_error("NaN detected in error measure");
        }
      m_error = *max_element(m_errors.begin(),m_errors.end());
      m_worst = max_element(m_errors.begin(),m_errors.end())-m_errors.begin();
      if (m_verbosity > 5) {
          xout << "IterativeSolverBase::calculateErrors m_errors"; for (size_t root=0; root<solution.size(); root++) xout <<" "<<m_errors[root]; xout <<std::endl;
        }
    }

    std::vector<double> IterativeSolverBase::eigenvalues()
    {
      std::vector<double> result;
      for (size_t root=0; root<(size_t)m_roots && root < (size_t)m_subspaceEigenvalues.rows(); root++) result.push_back(m_subspaceEigenvalues[root].real());
      return result;
    }

  };
}


#endif // ITERATIVESOLVER_H
