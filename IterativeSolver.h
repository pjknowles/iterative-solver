#ifndef ITERATIVESOLVER_H
#define ITERATIVESOLVER_H

#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <algorithm>
#include <stdexcept>
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
 template<class scalar=double>
 class IterativeSolver {
 public:
  IterativeSolver(
  ) :
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
    m_residual_eigen(false),
    m_residual_rhs(false),
    m_rspt(false),
    m_iterations(0),
    m_singularity_threshold(1e-20),
    m_error(0),
    m_lastVectorIndex(0),
    m_dimension(0),
    m_options(optionMap()) {}

  virtual ~IterativeSolver() = default;

  bool m_rspt;
 public:
  /*!
   * \brief Take, typically, a current solution and residual, and return new solution.
   * In the context of Lanczos-like linear methods, the input will be a current expansion vector and the result of
   * acting on it with the matrix, and the output will be a new expansion vector.
   * For non-linear equations, the input will be the current solution and residual, and the output the interpolated solution and residual.
   * \param parameters On input, the current solution or expansion vector. On exit, the interpolated solution vector.
   * \param action On input, the residual for parameters (non-linear), or action of matrix on parameters (linear). On exit, the expected (non-linear) or actual (linear) residual of the interpolated parameters.
   * \param parametersP On exit, the interpolated solution projected onto the P space.
   * \param other Optional additional vectors that should be interpolated like the residual.
   */
  void addVector(vectorSet<scalar> &parameters,
                 vectorSet<scalar> &action,
                 std::vector<std::vector<scalar> > &parametersP,
                 vectorSet<scalar> &other
  ) {
//   xout << "addVector initial parameters: "<<parameters<<std::endl;
   if (m_roots < 1) m_roots = parameters.size(); // number of roots defaults to size of parameters
   assert(parameters.size() == action.size());
   m_iterations++;
   for (size_t k = 0; k < action.size(); k++)
    action.m_active[k] = //action.m_active[k] &&
      parameters.m_active[k];
   m_lastVectorIndex = addVectorSet(parameters, action, other) -
                       1; // derivative classes might eventually store the vectors on top of previous ones, in which case they will need to store the position here for later calculation of iteration step
//   xout << "set lastVectorIndex=addVectorSet-1="<<m_lastVectorIndex<<std::endl;
   buildSubspace();
   solveReducedProblem();
   doInterpolation(parameters, action, parametersP, other);
  }

  void addVector(vectorSet<scalar> &parameters, vectorSet<scalar> &action) {
   std::vector<std::vector<scalar> > parametersP;
   return addVector(parameters, action, parametersP);
  }

  void
  addVector(vectorSet<scalar> &parameters, vectorSet<scalar> &action, std::vector<std::vector<scalar> > &parametersP) {
   vectorSet<scalar> other;
   return addVector(parameters, action, parametersP, other);
  }

 public:
  /*!
   * \brief Specify a P-space vector as a sparse combination of parameters. The container holds a number of segments,
   * each characterised by an offset in the full space, and a vector of coefficients starting at that offset.
   */
  using Pvector = std::map<size_t, scalar>;

  /*!
   * \brief Add P-space vectors to the expansion set for linear methods.
   * \param Pvectors the vectors to add
   * \param PP Matrix projected onto the existing+new, new P space. It should be provided as a
   * 1-dimensional array, with the existing+new index running fastest.
   * \param parameters On exit, the interpolated solution vector.
   * \param action On input, the residual for parameters (non-linear), or action of matrix on parameters (linear). On exit, the expected (non-linear) or actual (linear) residual of the interpolated parameters.
   * \param parametersP On exit, the interpolated solution projected onto the P space.
   * \param other On exit, interpolation of the other vectors
   */
  void addP(std::vector<Pvector> Pvectors, const scalar *PP,
            vectorSet<scalar> &parameters,
            vectorSet<scalar> &action,
            std::vector<std::vector<scalar> > &parametersP,
            vectorSet<scalar> &other) {
   size_t oldss = m_subspaceMatrix.rows();
   m_subspaceMatrix.conservativeResize(oldss + Pvectors.size(), oldss + Pvectors.size());
   m_subspaceOverlap.conservativeResize(oldss + Pvectors.size(), oldss + Pvectors.size());
   size_t old = m_PQMatrix.rows();
   parametersP.resize(old + Pvectors.size());
   m_PQMatrix.conservativeResize(old + Pvectors.size(), m_QQMatrix.rows());
   size_t offset = 0;
   for (size_t n = 0; n < Pvectors.size(); n++)
    m_Pvectors.push_back(Pvectors[n]);
   for (size_t n = 0; n < Pvectors.size(); n++) {
    for (int i = 0; i < m_subspaceMatrix.rows(); i++) {
     m_subspaceMatrix(old + n, i) = m_subspaceMatrix(i, old + n) = PP[offset++];
     double overlap = 0;
     for (const auto &p : Pvectors[n])
      if (m_Pvectors[i].count(p.first))
       overlap += p.second * m_Pvectors[i][p.first];
     m_subspaceOverlap(old + n, i) = m_subspaceOverlap(i, old + n) = overlap;
    }
   }
   size_t l = 0;
   for (size_t ll = 0; ll < m_solutions.size(); ll++) {
    for (size_t lll = 0; lll < m_solutions[ll].size(); lll++) {
     if (m_solutions[ll].m_active[lll]) {
      for (size_t n = 0; n < Pvectors.size(); n++) {
       m_PQMatrix(old + n, l) = m_residuals[ll][lll]->dot(Pvectors[n]);
       m_PQOverlap(old + n, l) = m_solutions[ll][lll]->dot(Pvectors[n]);
      }
      l++;
     }
    }
   }
   buildSubspace();
   solveReducedProblem();
   doInterpolation(parameters, action, parametersP, other);
  }

  void addP(std::vector<Pvector> Pvectors, const scalar *PP, vectorSet<scalar> &parameters, vectorSet<scalar> &action,
            std::vector<std::vector<scalar> > &parametersP) {
   vectorSet<scalar> other;
   return addP(Pvectors, PP, parameters, action, parametersP, other);
  }

  /*!
     * \brief Take the updated solution vector set, and adjust it if necessary so that it becomes the vector to
     * be used in the next iteration; this is done only in the case of linear solvers where the orthogonalize option is set.
     * Also calculate the degree of convergence, and write progress to std::cout.
     * \param solution The current solution, after interpolation and updating with the preconditioned residual.
     * \param residual The residual after interpolation.
     * \return true if convergence reached for all roots
     */
  bool endIteration(vectorSet<scalar> &solution, const vectorSet<scalar> &residual) {
   calculateErrors(solution, residual);
   if (m_error >= m_thresh) adjustUpdate(solution);
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
  std::vector<size_t> suggestP(const vectorSet<scalar> &solution,
                               const vectorSet<scalar> &residual,
                               const size_t maximumNumber = 1000,
                               const scalar threshold = 0) {
   std::map<size_t, scalar> result;
   for (size_t kkk = 0; kkk < solution.size(); kkk++) {
    if (solution.m_active[kkk]) {
     std::vector<size_t> indices;
     std::vector<scalar> values;
     std::tie(indices, values) = solution[kkk]->select(*residual[kkk], maximumNumber, threshold);
//     for (auto k=0; k<indices.size(); k++) std::cout << "select "<< indices[k] <<" : "<<values[k]<<std::endl;
     for (size_t i = 0; i < indices.size(); i++)
      if (result.count(indices[i]))
       result[indices[i]] = std::max(result[indices[i]], values[i]);
      else
       result[indices[i]] = values[i];
    }
   }
   // sort and select
//   for (const auto& kv : result) std::cout << "result: " << kv.first << " : " <<kv.second<<std::endl;
   std::multimap<scalar, size_t, std::greater<scalar> > inverseResult;
   for (const auto &kv : result) inverseResult.insert(std::pair<scalar, size_t>(kv.second, kv.first));
//   for (const auto& kv : inverseResult) std::cout << "inverseResult: " << kv.first << " : " <<kv.second<<std::endl;
   std::vector<size_t> indices;
//   std::vector<scalar> values;
   size_t k = 0;
   for (auto p = inverseResult.cbegin(); p != inverseResult.cend() && k < maximumNumber; k++) {
    indices.push_back(p->second);// values.push_back(p->first);
    ++p;
   }
//   for (auto k=0; k<indices.size(); k++) std::cout << "suggest P "<< indices[k] <<" : "<<values[k]<<std::endl;
   return indices;
  }


  /*!
   * \brief Set convergence threshold
   */
  void setThresholds(double thresh) { m_thresh = thresh; }

  unsigned int iterations() { return m_iterations; } //!< How many iterations have occurred

  std::vector<double> eigenvalues() ///< The calculated eigenvalues of m_subspaceMatrix
  {
   std::vector<double> result;
   for (size_t root = 0; root < (size_t) m_roots && root < (size_t) m_subspaceEigenvalues.rows(); root++)
    result.push_back(m_subspaceEigenvalues[root].real());
   return result;
  }

  std::vector<double> errors() const { return m_errors; } //!< Error at last iteration

  size_t dimensionP() const { return (size_t) m_PQMatrix.rows(); } //!< Size of P space

 private:
  Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> m_PQMatrix, m_PQOverlap; //!< The PQ block of the matrix
  std::vector<Pvector> m_Pvectors;
 public:
  int m_verbosity; //!< How much to print. Zero means nothing; One results in a single progress-report line printed each iteration.
  double m_thresh; //!< If predicted residual . solution is less than this, converged, irrespective of cthresh and gthresh.
  unsigned int m_maxIterations; //!< Maximum number of iterations in solve()
  unsigned int m_minIterations; //!< Minimum number of iterations in solve()
  bool m_orthogonalize; ///< Whether or not to orthogonalize the result of update() to all previous expansion vectors (appropriate only for linear methods).
  bool m_linear; ///< Whether residuals are linear functions of the corresponding expansion vectors.
  bool m_hermitian; ///< Whether residuals can be assumed to be the action of an underlying self-adjoint operator.
  size_t m_roots; ///< How many roots to calculate / equations to solve (defaults to size of solution and residual vectors)
 public:
  optionMap m_options; ///< A string of options to be interpreted by solveReducedProblem().

 private:
  void adjustUpdate(vectorSet<scalar> &solution) {
   //     xout << "m_errors[0] "<<m_errors[0]<<", m_thresh "<<m_thresh<<std::endl;
   for (size_t k = 0; k < solution.size(); k++)
    solution.m_active[k] = (m_errors[k] >= m_thresh || m_minIterations > m_iterations);
   if (m_orthogonalize) {
    //      xout << "IterativeSolverBase::adjustUpdate solution before orthogonalization: "<<solution<<std::endl;
    for (auto rep = 0; rep < 2; rep++)
     for (size_t kkk = 0; kkk < solution.size(); kkk++) {
      if (solution.m_active[kkk]) {
       for (auto i = 0; i < m_Pvectors.size(); i++) {
        const auto &p = m_Pvectors[i];
        double s = -solution[kkk]->dot(p) / m_subspaceOverlap(i, i);
        solution[kkk]->axpy(s, p);
       }
       for (size_t ll = 0; ll < m_solutions.size(); ll++) {
        for (size_t lll = 0; lll < m_solutions[ll].size(); lll++) {
         if (m_solutions[ll].m_active[lll]) {
          double s = -(m_solutions[ll][lll]->dot(*solution[kkk])) / (m_solutions[ll][lll]->dot(*m_solutions[ll][lll]));
          solution[kkk]->axpy(s, *m_solutions[ll][lll]);
         }
        }
       }
       for (size_t lll = 0; lll < kkk; lll++) {
        if (solution.m_active[lll]) {
         double s = solution[lll]->dot(*solution[kkk]);
         solution[kkk]->axpy(-s, *solution[lll]);
        }
       }
       double s = solution[kkk]->dot(*solution[kkk]);
       if (s <= 0)
        solution.m_active[kkk] = false;
       else
        solution[kkk]->scal(1 / std::sqrt(s));
      }
     }
    //      xout << "IterativeSolverBase::adjustUpdate solution after orthogonalization: "<<solution<<std::endl;
   }
  }

 protected:

  virtual void solveReducedProblem()=0;

  virtual void report() {
   if (m_verbosity > 0)
    xout << "iteration " << iterations() << ", error[" << m_worst << "] = " << m_error << std::endl;
  }

  void buildSubspace() {
   const auto nP = m_Pvectors.size();
   const auto nQ = m_QQMatrix.rows();
   const auto n = nP + nQ;
   m_subspaceMatrix.conservativeResize(n, n);
   m_subspaceOverlap.conservativeResize(n, n);
   for (auto i = 0; i < nQ; i++) {
    for (auto j = 0; j < nQ; j++) {
     m_subspaceMatrix(nP + j, nP + i) = m_QQMatrix(j, i);
     m_subspaceOverlap(nP + j, nP + i) = m_QQOverlap(j, i);
    }
    for (auto j = 0; j < nP; j++) {
     m_subspaceMatrix(j, nP + i) = m_subspaceMatrix(nP + i, j) = m_PQMatrix(j, i);
     m_subspaceOverlap(j, nP + i) = m_subspaceOverlap(nP + i, j) = m_PQOverlap(j, i);
    }
   }
   if (m_verbosity > 3 && m_PQMatrix.rows()>0) xout << "PQ matrix" << std::endl << this->m_PQMatrix << std::endl;
   if (m_verbosity > 3 && m_PQMatrix.rows()>0) xout << "PQ overlap" << std::endl << this->m_PQOverlap << std::endl;
   if (m_verbosity > 3 && m_PQMatrix.rows()>0) xout << "QQ matrix" << std::endl << this->m_QQMatrix << std::endl;
   if (m_verbosity > 3 && m_PQMatrix.rows()>0) xout << "QQ overlap" << std::endl << this->m_QQOverlap << std::endl;
   if (m_verbosity > 2) xout << "Subspace matrix" << std::endl << this->m_subspaceMatrix << std::endl;
   if (m_verbosity > 2) xout << "Subspace overlap" << std::endl << this->m_subspaceOverlap << std::endl;
  }

 protected:
  void diagonalizeSubspaceMatrix() {
   size_t kept = m_subspaceMatrix.rows();
   {
    Eigen::EigenSolver<Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> > ss(m_subspaceOverlap);
    Eigen::VectorXcd sse = ss.eigenvalues();
    for (int k = 0; k < sse.rows(); k++) {
     if (std::fabs(sse(k).real()) < m_singularity_threshold)
      kept--;
    }
   }
   if (m_verbosity >= 0 && kept < m_subspaceMatrix.rows())
    xout << "IterativeSolver WARNING, subspace singular, pruned from " << m_subspaceMatrix.rows() << " to " << kept
         << std::endl;

   Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> H = m_subspaceMatrix.block(0, 0, kept, kept);
   Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> S = m_subspaceOverlap.block(0, 0, kept, kept);
   Eigen::GeneralizedEigenSolver<Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic>> s(H, S);
   m_subspaceEigenvalues = s.eigenvalues();
   m_subspaceEigenvectors = s.eigenvectors();
   // sort
   std::vector<size_t> map;
   for (Eigen::Index k = 0; k < H.rows(); k++) {
    size_t ll;
    for (ll = 0; std::count(map.begin(), map.end(), ll) != 0; ll++);
    for (Eigen::Index l = 0; l < H.rows(); l++) {
     if (std::count(map.begin(), map.end(), l) == 0) {
      if (s.eigenvalues()(l).real() < s.eigenvalues()(ll).real())
       ll = l;
     }
    }
    map.push_back(ll);
    m_subspaceEigenvalues[k] = s.eigenvalues()(ll);
    for (Eigen::Index l = 0; l < H.rows(); l++) m_subspaceEigenvectors(l, k) = s.eigenvectors()(l, ll);
   }
   Eigen::MatrixXcd overlap = m_subspaceEigenvectors.transpose() * S * m_subspaceEigenvectors;
   for (Eigen::Index k = 0; k < overlap.rows(); k++)
    for (Eigen::Index l = 0; l < overlap.rows(); l++)
     m_subspaceEigenvectors(l, k) /= std::sqrt(overlap(k, k).real());
   //  xout << "eigenvalues"<<std::endl<<m_subspaceEigenvalues<<std::endl;
   //  xout << "eigenvectors"<<std::endl<<m_subspaceEigenvectors<<std::endl;
  }

  /*!
   * @brief form the combination of P and Q vectors to give the interpolated solution and corresponding residual (and maybe other vectors).
   * On entry, m_solution contains the interpolation
   *
   * @param solution On exit, the complete current solution (P and Q parts)
   * @param residual On exit, the Q contribution to the residual. The action of the matrix on the P solution is missing,
   * and has to be evaluated by the caller.
   * @param solutionP On exit, the solution projected to the P space
   * @param other On exit, interpolation of the other vectors
   */
  void doInterpolation(vectorSet<scalar> &solution, vectorSet<scalar> &residual,
                       std::vector<std::vector<scalar> > &solutionP,
                       vectorSet<scalar> &other) {
   solution.zero();
   residual.zero();
   size_t nP = m_Pvectors.size();
   solutionP.resize(residual.size());
   other.zero();
   for (size_t kkk = 0; kkk < residual.size() && kkk < m_interpolation.rows(); kkk++) {
    solutionP[kkk].resize(nP);
    for (size_t l = 0; l < nP; l++)
     solution[kkk]->axpy((solutionP[kkk][l] = this->m_interpolation(l, kkk)), m_Pvectors[l]);
    size_t l = nP;
    for (size_t ll = 0; ll < this->m_solutions.size(); ll++) {
     for (size_t lll = 0; lll < this->m_solutions[ll].size(); lll++) {
      if (this->m_solutions[ll].m_active[lll]) {
       if (m_verbosity > 2)
        xout << "LinearEigensystem::doInterpolation kkk=" << kkk << ", ll=" << ll << ", lll=" << lll << ", l=" << l
             << std::endl;
       if (m_verbosity > 2) xout << "Interpolation:\n" << this->m_interpolation(l, kkk) << std::endl;
       solution[kkk]->axpy(this->m_interpolation(l, kkk), *this->m_solutions[ll][lll]);
       residual[kkk]->axpy(this->m_interpolation(l, kkk), *this->m_residuals[ll][lll]);
       if (other.size() > kkk) other[kkk]->axpy(this->m_interpolation(l, kkk), *this->m_others[ll][lll]);
       l++;
      }
     }
    }
    if (m_residual_eigen) residual[kkk]->axpy(-this->m_subspaceEigenvalues(kkk).real(), *solution[kkk]);
    if (m_residual_rhs) residual[kkk]->axpy(-1, *this->m_rhs[kkk]);
   }

  }

  void calculateErrors(const vectorSet<scalar> &solution, const vectorSet<scalar> &residual) {
   if (m_verbosity > 5) {
    xout << "IterativeSolverBase::calculateErrors m_linear" << m_linear << std::endl;
    xout << "IterativeSolverBase::calculateErrors solution.m_active";
    for (size_t root = 0; root < solution.size(); root++) xout << " " << solution.m_active[root];
    xout << std::endl;
    xout << "IterativeSolverBase::calculateErrors solution " << solution << std::endl;
    xout << "IterativeSolverBase::calculateErrors residual.m_active";
    for (size_t root = 0; root < solution.size(); root++) xout << " " << residual.m_active[root];
    xout << std::endl;
    xout << "IterativeSolverBase::calculateErrors residual " << residual << std::endl;
   }
   m_errors.clear();
   if (m_linear) { // we can use the extrapolated residual if the problem is linear
    for (size_t k = 0; k < solution.size(); k++)
     m_errors.push_back(residual.m_active[k] ? std::fabs(residual[k]->dot(*solution[k])) : 0);
   } else {
    vectorSet<scalar> step = solution;
    step.axpy(-1, m_solutions[m_lastVectorIndex]);
    if (m_verbosity > 6) {
     xout << "IterativeSolverBase::calculateErrors last solution " << m_lastVectorIndex << std::endl;
     xout << "IterativeSolverBase::calculateErrors last solution " << m_solutions[m_lastVectorIndex] << std::endl;
     xout << "IterativeSolverBase::calculateErrors solution " << solution << std::endl;
     xout << "IterativeSolverBase::calculateErrors step " << step << std::endl;
    }
    for (size_t k = 0; k < solution.size(); k++)
     m_errors.push_back(
       m_residuals[m_lastVectorIndex].m_active[k] ? std::fabs(m_residuals[m_lastVectorIndex][k]->dot(*step[k])) : 1);
   }
   //  xout << "last active "<<m_lastVectorIndex<<" "<<m_residuals[m_lastVectorIndex].m_active[0]<<std::endl;
   for (const auto &e : m_errors) if (std::isnan(e)) throw std::overflow_error("NaN detected in error measure");
   m_error = *max_element(m_errors.begin(), m_errors.end());
   m_worst = max_element(m_errors.begin(), m_errors.end()) - m_errors.begin();
   if (m_verbosity > 5) {
    xout << "IterativeSolverBase::calculateErrors m_errors";
    for (size_t root = 0; root < solution.size(); root++) xout << " " << m_errors[root];
    xout << std::endl;
   }
  }

  size_t
  addVectorSet(const vectorSet<scalar> &solution, const vectorSet<scalar> &residual, const vectorSet<scalar> &other) {
   //      if (residual.m_active.front()==0) xout <<"warning: inactive residual"<<std::endl;
   //      if (solution.m_active.front()==0) xout <<"warning: inactive solution"<<std::endl;
   //      for (size_t kkk=0; kkk<solution.size(); kkk++)
   //          xout << "addVectorSet solution: "<<solution[kkk]<<std::endl;
   //      for (size_t kkk=0; kkk<residual.size(); kkk++)
   //          xout << "addVectorSet residual: "<<residual[kkk]<<std::endl;
   m_residuals.emplace_back(
     vectorSet<scalar>(residual, LINEARALGEBRA_CLONE_ADVISE_DISTRIBUTED | LINEARALGEBRA_CLONE_ADVISE_OFFLINE));
   m_solutions.emplace_back(
     vectorSet<scalar>(solution, LINEARALGEBRA_CLONE_ADVISE_DISTRIBUTED | LINEARALGEBRA_CLONE_ADVISE_OFFLINE));
   m_others.emplace_back(
     vectorSet<scalar>(other, LINEARALGEBRA_CLONE_ADVISE_DISTRIBUTED | LINEARALGEBRA_CLONE_ADVISE_OFFLINE));
   const vectorSet<scalar> &residual1 = residual;
   const vectorSet<scalar> &solution1 = solution;//      for (size_t kkk=0; kkk<solution.size(); kkk++)
   //          xout << "solution: "<<solution[kkk]<<std::endl;
   //      for (size_t kkk=0; kkk<residual.size(); kkk++)
   //          xout << "residual: "<<residual[kkk]<<std::endl;
   //      xout << "appendQQMatrix"<<std::endl;
   //      xout << "residual"<<std::endl<<residual[0]<<std::endl;
   //      xout << "solution"<<std::endl<<solution[0]<<std::endl;
   size_t nP = m_PQMatrix.rows();
   size_t old_size = m_QQMatrix.rows();
   size_t new_size = old_size + count(residual1.m_active.begin(), residual1.m_active.end(), true);
   //      xout << "old_size="<<old_size<<std::endl;
   //      xout << "new_size="<<new_size<<std::endl;
   m_QQMatrix.conservativeResize(new_size, new_size);
   m_QQOverlap.conservativeResize(new_size, new_size);
   m_PQMatrix.conservativeResize(nP, new_size);
   m_PQOverlap.conservativeResize(nP, new_size);
   std::vector<vectorSet<scalar>> *bra = m_subspaceMatrixResRes ? &m_residuals : &m_solutions;
   size_t k = old_size;
   for (size_t kkk = 0; kkk < residual1.size(); kkk++) {
    if (residual1.m_active[kkk]) {
     for (auto l = 0; l < nP; l++) {
      m_PQMatrix(l, k) = residual1[kkk]->dot(m_Pvectors[l]);
      m_PQOverlap(l, k) = solution1[kkk]->dot(m_Pvectors[l]);
     }
     size_t l = 0;
     for (size_t ll = 0; ll < m_solutions.size(); ll++) {
      for (size_t lll = 0; lll < m_solutions[ll].size(); lll++) {
       if (m_solutions[ll].m_active[lll]) {
//              xout << "bra"<<std::endl<<(*bra)[ll][lll]<<std::endl;
//              xout << "residual1"<<std::endl<<residual1[kkk]<<std::endl;
        m_QQMatrix(k, l) = m_QQMatrix(l, k) = (*bra)[ll][lll]->dot(*residual1[kkk]);
        m_QQOverlap(k, l) = m_QQOverlap(l, k) = m_solutions[ll][lll]->dot(*solution1[kkk]);
        l++;
       }
      }
     }
     m_dateOfBirth.push_back(++m_date);
     k++;
    }
   }
   if (m_verbosity > 3) {
    xout << "QQ Matrix: " << std::endl << m_QQMatrix << std::endl;
    xout << "QQ Overlap: " << std::endl << m_QQOverlap << std::endl;
   }
   return m_residuals.size();
  }

  void deleteVector(size_t index) {
   if (index >= m_residuals.size()) throw std::logic_error("invalid index");
   //    xout << "deleteVector "<<index<<std::endl;
   //    xout << "old m_subspaceMatrix"<<std::endl<<m_subspaceMatrix<<std::endl;
   //    xout << "old m_subspaceOverlap"<<std::endl<<m_subspaceOverlap<<std::endl;
   size_t old_size = m_QQMatrix.rows();
   size_t new_size = old_size;
   size_t l = 0;
   for (size_t ll = 0; ll < m_solutions.size(); ll++) {
    for (size_t lll = 0; lll < m_solutions[ll].size(); lll++) {
     if (m_solutions[ll].m_active[lll]) {
      if (l == index) { // this is the one to delete
       m_dateOfBirth.erase(m_dateOfBirth.begin() + l);
       m_solutions[ll].m_active[lll] = false;
       m_residuals[ll].m_active[lll] = false;
       //                    m_others[ll].m_active[lll] = false;
       for (int l2 = l + 1; l2 < m_QQMatrix.cols(); l2++) {
        for (int k = 0; k < m_QQMatrix.rows(); k++) {
         m_QQMatrix(k, l2 - 1) = m_QQMatrix(k, l2);
         m_QQOverlap(k, l2 - 1) = m_QQOverlap(k, l2);
        }
        for (int k = 0; k < m_PQMatrix.rows(); k++) {
         m_PQMatrix(k, l2 - 1) = m_PQMatrix(k, l2);
         m_PQOverlap(k, l2 - 1) = m_PQOverlap(k, l2);
        }
       }
       for (int l2 = l + 1; l2 < m_QQMatrix.rows(); l2++) {
        for (int k = 0; k < m_QQMatrix.rows(); k++) {
         m_QQMatrix(l2 - 1, k) = m_QQMatrix(l2, k);
         m_QQOverlap(l2 - 1, k) = m_QQOverlap(l2, k);
        }
       }
       new_size--;
      }
      l++;
     }
     if (std::count(m_solutions[ll].m_active.begin(), m_solutions[ll].m_active.end(), true) ==
         0) { // can delete the whole thing
     }
    }
   }
   m_QQMatrix.conservativeResize(new_size, new_size);
   m_QQOverlap.conservativeResize(new_size, new_size);
   //    xout << "new m_subspaceMatrix"<<std::endl<<m_subspaceMatrix<<std::endl;
   //    xout << "new m_subspaceOverlap"<<std::endl<<m_subspaceOverlap<<std::endl;
   buildSubspace();
  }

  std::vector<double> m_errors; //!< Error at last iteration
  double m_error; //!< worst error at last iteration
  size_t m_worst; //!< worst-converged solution, ie m_error = m_errors[m_worst]
  int m_date;
  bool m_subspaceMatrixResRes; // whether m_subspaceMatrix is Residual.Residual (true) or Solution.Residual (false)
  bool m_residual_eigen; // whether to subtract eigenvalue*solution when constructing residual
  bool m_residual_rhs; // whether to subtract rhs when constructing residual
  // whether to use RSPT to construct solution instead of diagonalisation
  std::vector<vectorSet<scalar>> m_residuals;
  std::vector<vectorSet<scalar>> m_solutions;
  std::vector<vectorSet<scalar>> m_others;
  vectorSet<scalar> m_rhs;
  std::vector<int> m_dateOfBirth;
  size_t m_lastVectorIndex;
  std::vector<scalar> m_updateShift;
  Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> m_interpolation; //!< The optimum combination of subspace vectors
  Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> m_subspaceMatrix;
  Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> m_subspaceOverlap;
  Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> m_QQMatrix;
  Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> m_QQOverlap;
  Eigen::MatrixXcd m_subspaceEigenvectors; // FIXME templating
  Eigen::VectorXcd m_subspaceEigenvalues; // FIXME templating
 public:
  size_t m_dimension;
 private:
  unsigned int m_iterations;
  double m_singularity_threshold;

 };

 template<class scalar>
 scalar inline
 operator*(const typename IterativeSolver<scalar>::Pvector &a, const typename IterativeSolver<scalar>::Pvector &b) {
  scalar result = 0;
  for (const auto &aa: a)
   if (b.find(aa.first))
    result += aa.second * b[aa.first];
 }
}

namespace LinearAlgebra {

 /*! @example LinearEigensystemExample.cpp */
 /*! @example LinearEigensystemExample-paged.cpp */
 /*!
 * \brief A class that finds the lowest eigensolutions of a matrix using Davidson's method, i.e. preconditioned Lanczos
 *
 * Example of simplest use with a simple in-memory container for eigenvectors: @include LinearEigensystemExample.cpp
 *
 * Example using a P-space and offline distributed storage provided by the PagedVector class: @include LinearEigensystemExample-paged.cpp
 *
 * \tparam scalar Type of matrix elements
 */
 template<class scalar=double>
 class LinearEigensystem : public IterativeSolver<scalar> {
 public:
  using IterativeSolver<scalar>::m_verbosity;

  /*!
   * \brief LinearEigensystem
   */
  LinearEigensystem() {
   this->m_residual_rhs = false;
   this->m_residual_eigen = true;
   this->m_linear = true;
   this->m_orthogonalize = true;
  }


 private:

  void solveReducedProblem() override {
   if (this->m_rspt) {
    throw std::logic_error("RSPT not yet implemented");
   } else {
    this->diagonalizeSubspaceMatrix();
    this->m_interpolation = this->m_subspaceEigenvectors.block(0, 0, this->m_subspaceEigenvectors.rows(),
                                                               std::min(int(this->m_roots), int(
                                                                 this->m_subspaceEigenvectors.rows()))).real();
   }

   this->m_updateShift.resize(this->m_roots);
   for (size_t root = 0; root < (size_t) this->m_roots; root++)
    this->m_updateShift[root] = -(1 + std::numeric_limits<scalar>::epsilon()) *
                                (root < this->m_subspaceEigenvectors.rows() ? this->m_subspaceEigenvalues[root].real()
                                                                            : 0);
  }

  void report() override {
   std::vector<scalar> ev = this->eigenvalues();
   if (m_verbosity > 0) {
    xout << "iteration " << this->iterations() << ", error[" << this->m_worst << "] = " << this->m_error
         << ", eigenvalues: ";
    for (const auto e : ev) xout << " " << e;
    xout << std::endl;
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
 template<class scalar=double>
 class LinearEquations : public IterativeSolver<scalar> {
 public:
  using IterativeSolver<scalar>::m_verbosity;

  /*!
   * \brief Constructor
   * \param rhs right-hand-side vectors. More can be added subsequently using addEquations(), provided iterations have not yet started.
   */
  LinearEquations(const vectorSet<scalar> &rhs)
    : IterativeSolver<scalar>::m_linear(true), IterativeSolver<scalar>::m_orthogonalize(true) {
   this->m_residual_rhs = true;
   this->m_residual_eigen = false;
   addEquations(rhs);
  }

  LinearEquations()
    : IterativeSolver<scalar>::m_linear(true), IterativeSolver<scalar>::m_orthogonalize(true) {}

  /*!
   * \brief add one or more equations to the set to be solved, by specifying their right-hand-side vector
   * \param rhs right-hand-side vectors to be added
   */
  void addEquations(const vectorSet<scalar> &rhs) {
   for (const auto &v : rhs) this->m_rhs.push_back(v);
  }


 protected:
  virtual void solveReducedProblem() {
   throw std::runtime_error("Not yet coded");
   this->diagonalizeSubspaceMatrix();
   this->m_interpolation.setZero();
  }

  virtual void report() {
   throw std::runtime_error("Not yet coded");
   std::vector<scalar> ev = this->eigenvalues();
   if (m_verbosity > 0) {
    xout << "iteration " << this->iterations() << ", error[" << this->m_worst << "] = " << this->m_error
         << ", eigenvalues: ";
    for (const auto e : ev) xout << " " << e;
    xout << std::endl;
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
 template<class scalar=double>
 class DIIS : public IterativeSolver<scalar> {
  using IterativeSolver<scalar>::m_residuals;
  using IterativeSolver<scalar>::m_solutions;
  using IterativeSolver<scalar>::m_others;
 public:
  using IterativeSolver<scalar>::m_verbosity;
  enum DIISmode_type {
   disabled ///< No extrapolation is performed
   , DIISmode ///< Direct Inversion in the Iterative Subspace
   , KAINmode ///< Krylov Accelerated Inexact Newton
  };

  /*!
 * \brief DIIS
 */
  DIIS()
    : m_svdThreshold(1e-10), m_maxDim(6) {
   this->m_residual_rhs = false;
   this->m_residual_eigen = false;
   this->m_orthogonalize = false;
   this->m_roots = 1;
   setMode(DIISmode);
   Reset();
  }

  ~DIIS() {}


  /*!
 * \brief Set options for DIIS.
 * \param mode Whether to perform DIIS, KAIN, or nothing.
 */
  virtual void setMode(enum DIISmode_type mode = DIISmode) {
   m_DIISmode = mode;
   this->m_subspaceMatrixResRes = mode != KAINmode;
//     this->m_preconditionResiduals = mode==KAINmode; // FIXME

   Reset();
   if (m_verbosity > 1)
    xout << "m_DIISmode set to " << m_DIISmode << std::endl;
  }

  /*!
 * \brief discards previous iteration vectors, but does not clear records
 */
  void Reset() {
   m_LastResidualNormSq = 1e99; // so it can be tested even before extrapolation is done
   this->m_lastVectorIndex = 0;
   while (this->m_QQMatrix.rows() > 0) this->deleteVector(0);;
   this->m_residuals.clear();
   this->m_solutions.clear();
   this->m_others.clear();
   this->m_Weights.clear();
  }
 protected:
  void solveReducedProblem() override {
   //	  xout << "Enter DIIS::solveReducedProblem"<<std::endl;
   //	  xout << "residual : "<<residual<<std::endl;
   //	  xout << "solution : "<<solution<<std::endl;
   this->m_updateShift.clear();
   this->m_updateShift.push_back(-(1 + std::numeric_limits<double>::epsilon()) * this->m_subspaceMatrix(0, 0));
   double weight =
     this->m_options.count("weight") ? (
       this->m_options.find("weight")->second) : 1.0;
   if (this->m_maxDim <= 1 || this->m_DIISmode == disabled) return;

   if (this->m_roots > 1) throw std::logic_error("DIIS does not handle multiple solutions");

   //  if (m_subspaceMatrix.rows() < 9) {
   //      xout << "m_subspaceMatrix on entry to DIIS::solveReducedProblem"<<std::endl<<m_subspaceMatrix<<std::endl;
   //  }
   size_t nDim = this->m_subspaceMatrix.rows();
   this->m_LastResidualNormSq = std::fabs(this->m_subspaceMatrix(nDim - 1, nDim - 1));
   //  xout << "this->m_LastResidualNormSq "<<this->m_LastResidualNormSq<<std::endl;

   this->m_Weights.push_back(weight); // TODO not conformant - add style

   Eigen::Array<scalar, Eigen::Dynamic, Eigen::Dynamic> d = this->m_subspaceMatrix.diagonal().array().abs();
   int worst = 0, best = 0;
   for (size_t i = 0; i < nDim; i++) {
    if (d(i) > d(worst)) worst = i;
    if (d(i) < d(best)) best = i;
   }
   double fBaseScale = std::sqrt(d(worst) * d(best));

   while (nDim > this->m_maxDim) { // prune away the worst/oldest vector. Algorithm to be done properly yet
    size_t prune = worst;
    if (true || prune == nDim - 1) { // the current vector is the worst, so delete the oldest
     //          xout << "this->m_dateOfBirth: "; for (auto b=this->m_dateOfBirth.begin(); b!=this->m_dateOfBirth.end(); b++) xout <<(*b); xout<<std::endl;
     prune = std::min_element(this->m_dateOfBirth.begin(), this->m_dateOfBirth.end()) - this->m_dateOfBirth.begin();
    }
    //      xout << "prune="<<prune<<std::endl;
    //  xout << "nDim="<<nDim<<", m_Weights.size()="<<m_Weights.size()<<std::endl;
    this->deleteVector(prune);
    m_Weights.erase(m_Weights.begin() + prune);
    nDim--;
    //  xout << "nDim="<<nDim<<", m_Weights.size()="<<m_Weights.size()<<std::endl;
   }
   if (nDim != (size_t) this->m_subspaceMatrix.rows()) throw std::logic_error("problem in pruning");
   if (m_Weights.size() != (size_t) this->m_subspaceMatrix.rows()) {
    xout << "nDim=" << nDim << ", m_Weights.size()=" << m_Weights.size() << std::endl;
    throw std::logic_error("problem after pruning weights");
   }



   //  if (m_subspaceMatrix.rows() < 9) {
   //      xout << "m_subspaceMatrix"<<std::endl<<m_subspaceMatrix<<std::endl;
   //  }

   // build actual DIIS system for the subspace used.
   Eigen::VectorXd
     Rhs(nDim + 1),
     Coeffs(nDim + 1);
   Eigen::MatrixXd
     B(nDim + 1, nDim + 1);

   // Factor out common size scales from the residual dots.
   // This is done to increase numerical stability for the case when _all_
   // residuals are very small.
   B.block(0, 0, nDim, nDim) = this->m_subspaceMatrix / fBaseScale;
   Rhs.head(nDim) = Eigen::VectorXd::Zero(nDim);

   // make Lagrange/constraint lines.
   for (size_t i = 0; i < nDim; ++i)
    B(i, nDim) = B(nDim, i) = -m_Weights[i];
   B(nDim, nDim) = 0.0;
   Rhs[nDim] = -1;
   //  xout << "B:"<<std::endl<<B<<std::endl;
   //  xout << "Rhs:"<<std::endl<<Rhs<<std::endl;

   // invert the system, determine extrapolation coefficients.
   Eigen::JacobiSVD<Eigen::MatrixXd> svd(B, Eigen::ComputeThinU | Eigen::ComputeThinV);
   svd.setThreshold(m_svdThreshold);
   Coeffs = svd.solve(Rhs).head(nDim);
   m_LastAmplitudeCoeff = Coeffs[nDim - 1];
   if (m_verbosity > 1) xout << "Combination of iteration vectors: " << Coeffs.transpose() << std::endl;
   for (size_t k = 0; k < (size_t) Coeffs.rows(); k++)
    if (std::isnan(Coeffs(k))) {
     xout << "B:" << std::endl << B << std::endl;
     xout << "Rhs:" << std::endl << Rhs << std::endl;
     xout << "Combination of iteration vectors: " << Coeffs.transpose() << std::endl;
     throw std::overflow_error("NaN detected in DIIS submatrix solution");
    }
   this->m_interpolation = Coeffs.head(nDim);
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

  unsigned int nLastDim() const { return m_residuals.size(); }

  unsigned int nMaxDim() const { return m_maxDim; }

  double m_svdThreshold; ///< Threshold for singular-value truncation in linear equation solver.
  size_t m_maxDim; ///< Maximum DIIS dimension allowed.
  static void
  randomTest(size_t sample, size_t n = 100, double alpha = 0.1, double gamma = 0.0, DIISmode_type mode = DIISmode);

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
}

// C interface
extern "C" void
IterativeSolverLinearEigensystemInitialize(size_t nQ, size_t nroot, double thresh, unsigned int maxIterations, int verbosity);

extern "C" void
IterativeSolverDIISInitialize(size_t n, double thresh, unsigned int maxIterations, int verbosity);

extern "C" void IterativeSolverFinalize();

extern "C" void
IterativeSolverAddVector(double *parameters, double *action, double *parametersP);

extern "C" int IterativeSolverEndIteration(double *c, double *g, double *error);

extern "C" void IterativeSolverAddP(const size_t nP, const size_t *offsets, const size_t *indices,
                                    const double *coefficients, const double *pp,
                                    double *parameters, double *action, double *parametersP);

extern "C" void IterativeSolverEigenvalues(double *eigenvalues);


#endif // ITERATIVESOLVER_H
