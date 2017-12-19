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

  typedef double scalar;
  typedef vector<scalar> ParameterVector;
  typedef vectorSet<scalar> ParameterVectorSet;
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
  class IterativeSolverBase
  {
  public:
    //    typedef void (*ParameterSetTransformation)(const ParameterVectorSet & inputs, ParameterVectorSet & outputs, std::vector<scalar> shift, bool append);
    struct ParameterSetTransformation
    {
      virtual void operator()(const ParameterVectorSet & inputs, ParameterVectorSet & outputs, std::vector<LinearAlgebra::scalar> shift, bool append) const =0;
    };
    /*!
 * \brief Place-holding template for residual calculation. It just returns the input as output.
 * \param inputs The parameters.
 * \param outputs On output, contains the corresponding residual vectors.
 * \param shift
 * \param append Whether to add the result to the original content of outputs
 */
    static struct : ParameterSetTransformation {
      void operator()(const ParameterVectorSet & inputs, ParameterVectorSet & outputs, std::vector<scalar> shift=std::vector<scalar>(), bool append=false) const override { if (append) outputs=inputs; }
    } noOp ;
    /*!
 * \brief Place-holding template for update calculation. It just returns the input as output.
 * \param inputs The parameters.
 * \param append Whether to add the result to the original content of outputs
 * \param outputs On output, contains the corresponding residual vectors.
 * \param shift
 */
    static struct : ParameterSetTransformation {
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
    IterativeSolverBase(const ParameterSetTransformation& residualFunction, const ParameterSetTransformation& preconditionerFunction=steepestDescent);
    virtual ~IterativeSolverBase();
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
    virtual bool iterate(ParameterVectorSet & residual, ParameterVectorSet & solution, ParameterVectorSet & other, const optionMap options=optionMap());
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
  };
}


#endif // ITERATIVESOLVER_H
