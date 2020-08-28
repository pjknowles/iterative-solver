#ifndef ITERATIVESOLVER_H
#define ITERATIVESOLVER_H
#include <vector>
#include <iostream>
#include "molpro/linalg/iterativesolver/P.h"
#include "molpro/linalg/iterativesolver/Q.h"
#include "molpro/linalg/iterativesolver/Statistics.h"
#include "molpro/linalg/iterativesolver/helper.h"
#include <molpro/linalg/iterativesolver/ArrayHandlers.h>
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif

#include <algorithm>
#include <cassert>
#include <climits>
#include <cmath>
#include <complex>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifndef LINEARALGEBRA_OFFLINE
#define LINEARALGEBRA_OFFLINE 0x01
#endif
#ifndef LINEARALGEBRA_DISTRIBUTED
#define LINEARALGEBRA_DISTRIBUTED 0x02
#endif

#undef isnan
#undef isinf
#include <memory>
#include <molpro/iostream.h>

/*!
 * @brief Contains classes that implement various iterative equation solvers.
 * They all share the feature that access to the client-provided solution and residual
 * vectors is via a potentially opaque interface to copy, scale, scalar product and
 * scalar-times-vector operations only.
 */
namespace molpro {
namespace linalg {
typedef std::map<std::string, std::string> optionMap;
template <class T>
static std::vector<typename T::value_type> nullVectorP;
template <class T>
static std::vector<std::vector<typename T::value_type>> nullVectorSetP;
template <class T>
static std::vector<std::reference_wrapper<std::vector<typename T::value_type>>> nullVectorRefSetP;

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
 * appropriate preconditioning on the residual, if addVector() has requested it.
 * -  optionally (but necessarily for BFGS optimization) make a call to endIteration(). Except in the BFGS case,
 * endIteration() merely calls report(), which prints progress information to standard output
 *
 * Classes that derive from this will, in the simplest case, need to provide just the solveReducedProblem() method that
 * governs how the parameter and residual vectors from successive iterations should be combined to form an optimum
 * solution with minimal residual.
 *
 * The underlying vector spaces are accessed through instances of the Rvector class.
 * - Both scalar (Rvector) and vector(std::vector<Rvector>) interfaces are provided to the class functions addVector(),
 * addValue() and endIteration()
 * - The following operations are carried out on Rvector objects
 *   - scalar_type dot(const Rvector& x, const Rvector& y) - scalar product of two vectors
 *   - void axpy(scalar_type a, Rvector& x, const Rvector& y) add a multiple a of y to x;
 *   - scal(Rvector& x, scalar_type a) scale x by a. In the special case a==0, x on entry might be uninitialised.
 * and these operations need to be provided by a handler class.
 *
 * @tparam Rvector The class encapsulating solution and residual vectors passed through the interface to class
 * functions.
 * @tparam Qvector Used internally as a class for storing vectors on backing store.
 * Handler classes must be provided for operations between Rvector and Qvector classes. Qvector objects are created by
 * copy-construction from an Rvector object, and no functions that subsequently modify them are called.
 * @tparam Pvector Specify a P-space vector as a sparse combination of parameters.
 */
template <class Rvector = std::vector<double>, class Qvector = Rvector,
          class Pvector = std::map<size_t, typename Rvector::value_type>>
class IterativeSolver {
public:
  IterativeSolver(const std::shared_ptr<iterativesolver::ArrayHandlers<Rvector, Qvector, Pvector>>& handlers)
      : // clang-format off
      m_handlers(handlers),
      m_verbosity(0),
      m_thresh(1e-8),
      m_maxIterations(1000),
      m_minIterations(0),
      m_linear(false),
      m_hermitian(false),
      m_roots(0),
      m_rspt(false),
      m_options(optionMap()),
      m_subspaceMatrixResRes(false),
      m_residual_eigen(false),
      m_residual_rhs(false),
      m_rhs(),
      m_lastVectorIndex(0),
      m_updateShift(0),
      m_dimension(0),
      m_value_print_name("value"),
      m_singularity_threshold(1e-4),
      m_augmented_hessian(0),
      m_svdThreshold(1e-15),
      m_maxQ(std::max(m_roots, size_t(16))),
      m_pspace(),
      m_qspace(m_pspace, m_hermitian, m_handlers),
      m_exclude_r_from_redundancy_test(false),
      m_orthogonalise_Q(true),
      m_nullify_solution_before_update(false)
{}
  // clang-format on

  IterativeSolver() = delete;
  IterativeSolver(const IterativeSolver&) = default;
  IterativeSolver(IterativeSolver&&) noexcept = default;
  IterativeSolver& operator=(const IterativeSolver&) = default;
  IterativeSolver& operator=(IterativeSolver&&) noexcept = default;
  virtual ~IterativeSolver() = default;

protected:
  using value_type = typename Rvector::value_type; ///< The underlying type of elements of vectors
  using vectorSet = typename std::vector<Rvector>; ///< Container of vectors
  using vectorRefSet = typename std::vector<std::reference_wrapper<Rvector>>;            ///< Container of vectors
  using constVectorRefSet = typename std::vector<std::reference_wrapper<const Rvector>>; ///< Container of vectors
  using vectorP = typename std::vector<value_type>;                                      ///< P-space parameters
  using vectorRefSetP = typename std::vector<std::reference_wrapper<vectorP>>; ///< Container of P-space parameters
  using constVectorRefSetP =
      typename std::vector<std::reference_wrapper<const vectorP>>; ///< Container of P-space parameters
  using vectorSetP = typename std::vector<vectorP>;                ///< Container of P-space parameters
  mutable std::shared_ptr<iterativesolver::ArrayHandlers<Rvector, Qvector, Pvector>> m_handlers;

public:
  std::shared_ptr<iterativesolver::ArrayHandlers<Rvector, Qvector, Pvector>> handlers() { return m_handlers; };
  using scalar_type =
      typename array::ArrayHandler<Rvector, Qvector>::value_type; ///< The type of scalar products of vectors
  /*!
   * \brief Take, typically, a current solution and residual, and return new solution.
   * In the context of Lanczos-like linear methods, the input will be a current expansion vector and the result of
   * acting on it with the matrix, and the output will be a new expansion vector.
   * For non-linear equations, the input will be the current solution and residual, and the output the interpolated
   * solution and residual. \param parameters On input, the current solution or expansion vector. On exit, the
   * interpolated solution vector. \param action On input, the residual for parameters (non-linear), or action of matrix
   * on parameters (linear). On exit, the expected (non-linear) or actual (linear) residual of the interpolated
   * parameters. \param parametersP On exit, the interpolated solution projected onto the P space.
   * \return whether it is expected that the client should make an update, based on the returned parameters and
   * residual, before the subsequent call to endIteration()
   */
  size_t addVector(vectorRefSet parameters, vectorRefSet action,
                   vectorRefSetP parametersP = nullVectorRefSetP<Rvector>) {
    if (m_roots < 1)
      m_roots = parameters.size();                      // number of roots defaults to size of parameters
    if (m_qspace.size() == 0 and m_working_set.empty()) // initial
      for (auto i = 0; i < parameters.size(); i++)
        m_working_set.push_back(i);
    //    molpro::cout << "m_working_set size " << m_working_set.size() << std::endl;
    if (m_working_set.size() == 0)
      return 0;
    assert(parameters.size() >= m_working_set.size());
    assert(parameters.size() == action.size());
    m_statistics.iterations++;
    m_statistics.r_creations += m_working_set.size();
    m_current_r.clear();
    m_current_v.clear();
    for (size_t k = 0; k < m_working_set.size(); k++) {
      if (m_residual_eigen) { // scale to roughly unit length for homogeneous equations in case the update has produced
                              // a very large vector in response to degeneracy
        auto s = m_handlers->rr().dot(parameters[k], parameters[k]);
        if (std::abs(s - 1) > 1e-3) {
          m_handlers->rr().scal(1 / std::sqrt(s), parameters[k]);
          m_handlers->rr().scal(1 / std::sqrt(s), action[k]);
        }
      }
      m_statistics.current_r_creations++;
      m_current_r.emplace_back(m_handlers->qr().copy(parameters[k]));
      m_current_v.emplace_back(m_handlers->qr().copy(action[k]));
    }
    if (not m_last_d.empty()) {
      assert(m_last_d.size() == m_working_set.size());
      assert(m_last_hd.size() == m_working_set.size());
      for (size_t k = 0; k < m_working_set.size(); k++) {
        m_statistics.q_creations++;
        m_qspace.add(parameters[k], action[k], m_last_d[k], m_last_hd[k], m_rhs, m_subspaceMatrixResRes,
                     m_orthogonalise_Q);
      }
      m_last_d.clear();
      m_last_hd.clear();
    }
    // TODO this generates another read for the q space which could perhaps be avoided
    for (size_t a = 0; a < m_qspace.size(); a++) {
      m_s_qr[a] = std::vector<value_type>(m_working_set.size());
      m_h_qr[a] = std::vector<value_type>(m_working_set.size());
      m_hh_qr[a] = std::vector<value_type>(m_working_set.size());
      m_h_rq[a] = std::vector<value_type>(m_working_set.size());
      const auto& qa = m_qspace[a];
      const auto& qha = m_qspace.action(a);
      for (size_t m = 0; m < m_working_set.size(); m++) {
        m_s_qr[a][m] = m_handlers->rq().dot(parameters[m], m_qspace[a]);
        m_h_qr[a][m] = m_handlers->rq().dot(action[m], m_subspaceMatrixResRes ? m_qspace.action(a) : m_qspace[a]);
        m_hh_qr[a][m] = m_handlers->rq().dot(action[m], m_qspace.action(a));
        m_h_rq[a][m] = m_hermitian ? m_h_qr[a][m]
                                   : m_handlers->rq().dot((m_subspaceMatrixResRes ? action[m] : parameters[m]),
                                                          m_qspace.action(a));
        //        molpro::cout << "a=" << a << ", m=" << m << ", m_s_qr " << m_s_qr[a][m] << ", m_h_qr " << m_h_qr[a][m]
        //                     << ", m_h_rq " << m_h_rq[a][m] << std::endl;
      }
    }
    m_s_pr.clear();
    m_h_pr.clear();
    m_h_rp.clear();
    for (auto p = 0; p < m_pspace.size(); p++) {
      m_s_pr.push_back(std::vector<value_type>(m_working_set.size()));
      m_h_pr.push_back(std::vector<value_type>(m_working_set.size()));
      m_h_rp.push_back(std::vector<value_type>(m_working_set.size()));
      for (size_t k = 0; k < m_working_set.size(); k++) {
        m_s_pr[p][k] = m_handlers->rp().dot(parameters[k].get(), m_pspace[p]);
        m_h_pr[p][k] = m_h_rp[p][k] =
            m_handlers->rp().dot(action[k].get(),
                                 m_pspace[p]); // TODO this works only for hermitian. Check that there is a
        // check
      }
    }
    m_s_rr.clear();
    m_h_rr.clear();
    m_hh_rr.clear();
    m_rhs_r.clear();
    for (auto m = 0; m < m_working_set.size(); m++) {
      m_rhs_r.push_back(std::vector<value_type>(m_rhs.size()));
      m_s_rr.push_back(std::vector<value_type>(m_working_set.size()));
      m_h_rr.push_back(std::vector<value_type>(m_working_set.size()));
      m_hh_rr.push_back(std::vector<value_type>(m_working_set.size()));
      for (size_t rhs = 0; rhs < m_rhs.size(); rhs++)
        m_rhs_r[m][rhs] = m_handlers->rr().dot(parameters[m], m_rhs[rhs]);
      for (size_t n = 0; n < m_working_set.size(); n++) {
        m_s_rr[m][n] = m_handlers->rr().dot(parameters[n], parameters[m]);
        m_h_rr[m][n] = m_handlers->rr().dot(action[n], (m_subspaceMatrixResRes ? action[m] : parameters[m]));
        m_hh_rr[m][n] = m_handlers->rr().dot(action[n], action[m]);
      }
    }

    return solveAndGenerateWorkingSet(parameters, action, parametersP);
  }
  size_t solveAndGenerateWorkingSet(vectorRefSet parameters, vectorRefSet action,
                                    vectorRefSetP parametersP = nullVectorRefSetP<Rvector>,
                                    bool calculateError = true) {
    buildSubspace();
    solveReducedProblem();
    //    molpro::cout << "update=" << update << std::endl;
    //    calculateResidualConvergence();
    // move any newly-converged solutions into the q space
    //      molpro::cout << "errors after calculateResidualConvergence()";
    //      for (const auto& e : m_errors)
    //        molpro::cout << " " << e;
    //      molpro::cout << std::endl;
    // evaluate all residuals
    // TODO make this more efficient by doing action just once, and not, maybe, for converged roots
    if (m_roots > parameters.size()) // TODO remove this restriction
      throw std::runtime_error("Cannot yet work with buffer smaller than number of roots");
    m_errors.resize(m_roots);
    m_working_set.clear();
    for (auto root = 0; root < m_roots; root++)
      m_working_set.push_back(root);
    if (calculateError) {
      if (m_linear)
        doInterpolation(parameters, action, parametersP, false);
      for (auto k = 0; k < m_working_set.size(); k++) {
        double a = std::abs(m_handlers->rr().dot(action[k], action[k]));
        m_errors[m_working_set[k]] = std::sqrt(a);
      }
    } else
      m_errors.assign(m_roots, 1); // TODO not right for retired roots

    doInterpolation(parameters, action, parametersP, true);
    m_last_d.clear();
    m_last_hd.clear();
    //    molpro::cout << "working set size "<<m_working_set.size()<<std::endl;
    //    molpro::cout << "m_thresh "<<m_thresh<<std::endl;
    for (int k = 0; k < m_working_set.size(); k++) {
      auto root = m_working_set[k];
      //      molpro::cout << "k=" << k << ", root=" << root << ", error=" << m_errors[root] << std::endl;
      if (m_linear and m_errors[root] < m_thresh and m_q_solutions.count(root) == 0) { // converged just now
        if (m_verbosity > 1)
          molpro::cout << "selecting root " << root << " for adding converged solution to Q space at position "
                       << m_qspace.size() << std::endl;
        m_statistics.q_creations++;
        m_qspace.add(parameters[k], action[k], m_rhs, m_subspaceMatrixResRes);
        m_q_solutions[m_working_set[k]] = m_qspace.keys().back();
      }
      if (m_linear and m_errors[root] < m_thresh) { // converged
        //        molpro::cout << "  remove this vector from the working set"<<std::endl;
        //  remove this vector from the working set
        // FIXME Doesn't this cause a copy? Should use swap, without .get
        for (auto kp = k + 1; kp < m_working_set.size(); kp++) {
          parameters[kp - 1].get() = parameters[kp].get();
          action[kp - 1].get() = action[kp].get();
          m_working_set[kp - 1] = m_working_set[kp];
        }
        m_working_set.pop_back();
        k--;
        //        molpro::cout << "k now = "<<k<<std::endl;
        //        molpro::cout << "m_working_set revised to ";for (const auto& w : m_working_set) molpro::cout <<"
        //        "<<w;molpro::cout <<std::endl;
      } else { // unconverged
               //        molpro::cout << "unconverged"<<std::endl;
        m_statistics.d_creations++;
        m_last_d.emplace_back(m_handlers->qr().copy(parameters[k]));
        m_last_hd.emplace_back(m_handlers->qr().copy(action[k]));
      }
    }
    //    molpro::cout << "working set size "<<m_working_set.size()<<std::endl;
    //    molpro::cout << "m_last_d size "<<m_last_d.size()<<std::endl;
    assert(m_last_d.size() == m_working_set.size());

    // re-establish the residual
    // TODO make more efficient
    doInterpolation(parameters, action, parametersP, false);
    if (m_nullify_solution_before_update) {
      m_last_d.clear();
      m_last_hd.clear();
      for (auto k = 0; k < m_working_set.size(); k++) {
        m_handlers->rr().fill(0, parameters[k].get());
        // FIXME Does this require a copy or is it a move?
        m_statistics.d_creations++;
        m_last_d.emplace_back(m_current_r[k]);
        m_last_hd.emplace_back(m_current_v[k]);
      }
    }
    m_current_r.clear();
    m_current_v.clear();
    return m_working_set.size();
  }
  size_t addVector(std::vector<Rvector>& parameters, std::vector<Rvector>& action,
                   vectorSetP& parametersP = nullVectorSetP<Rvector>) {
    return addVector(vectorRefSet(parameters.begin(), parameters.end()), vectorRefSet(action.begin(), action.end()),
                     vectorRefSetP(parametersP.begin(), parametersP.end()));
  }
  size_t addVector(Rvector& parameters, Rvector& action, vectorP& parametersP = nullVectorP<Rvector>) {
    return addVector(vectorRefSet(1, parameters), vectorRefSet(1, action), vectorRefSetP(1, parametersP));
  }

  /*!
   * \brief Take a current solution, objective function value and residual, and return new solution.
   * \param parameters On input, the current solution. On exit, the interpolated solution vector.
   * \param value The value of the objective function for parameters.
   * \param action On input, the residual for parameters. On exit, the expected (non-linear) residual of the
   * interpolated parameters. \return whether it is expected that the client should make an update, based on the
   * returned parameters and residual, before the subsequent call to endIteration()
   */
  size_t addValue(Rvector& parameters, value_type value, Rvector& action) {
    m_values.push_back(value);
    return this->addVector(parameters, action);
  }

public:
  /*!
   * \brief Add P-space vectors to the expansion set for linear methods.
   * \param Pvectors the vectors to add
   * \param PP Matrix projected onto the existing+new, new P space. It should be provided as a
   * 1-dimensional array, with the existing+new index running fastest.
   * \param parameters On exit, the interpolated solution vector.
   * \param action  On exit, the  residual of the interpolated Q parameters.
   * The contribution from the new, and any existing, P parameters is missing, and should be added in subsequently.
   * \param parametersP On exit, the interpolated solution projected onto the P space.
   * \return The number of vectors contained in parameters, action, parametersP
   */
  size_t addP(std::vector<Pvector> Pvectors, const value_type* PP, vectorRefSet parameters, vectorRefSet action,
              vectorRefSetP parametersP) {
    m_statistics.p_creations += Pvectors.size();
    m_pspace.add(Pvectors, PP, m_rhs, m_handlers->pp(), m_handlers->qp());
    m_qspace.refreshP(action.front());
    m_working_set.clear();
    auto result = solveAndGenerateWorkingSet(parameters, action, parametersP, false);
    m_last_d.clear(); // TODO more intelligent way needed
    m_last_hd.clear();
    return result;
  }

  size_t addP(std::vector<Pvector> Pvectors, const value_type* PP, std::vector<Rvector>& parameters,
              std::vector<Rvector>& action, vectorSetP& parametersP) {
    return addP(Pvectors, PP, vectorRefSet(parameters.begin(), parameters.end()),
                vectorRefSet(action.begin(), action.end()), vectorRefSetP(parametersP.begin(), parametersP.end()));
  }
  size_t addP(Pvector Pvectors, const value_type* PP, Rvector& parameters, Rvector& action, vectorP& parametersP) {
    return addP(std::vector<Pvector>{Pvectors}, PP, vectorRefSet(1, parameters), vectorRefSet(1, action),
                vectorRefSetP(1, parametersP));
  }

  /*!
   * @brief Empty the P space
   */
  void clearP() { m_qspace.clearP(); }

  /*!
   * \brief For most solvers, this function does nothing but report, but the exception is Optimize.
   * Also write progress to molpro::cout.
   * \param solution The current
   * solution, after interpolation and updating with the preconditioned residual.
   * \param residual The residual after interpolation.
   * \return Whether convergence has been reached
   */
  virtual bool endIteration(vectorRefSet solution, constVectorRefSet residual) {
    if (m_verbosity >= 0)
      report();
    return m_errors.empty() ? false : *std::max_element(m_errors.cbegin(), m_errors.cend()) < m_thresh;
  }
  virtual bool endIteration(std::vector<Rvector>& solution, const std::vector<Rvector>& residual) {
    return endIteration(vectorRefSet(solution.begin(), solution.end()),
                        constVectorRefSet(residual.begin(), residual.end()));
  }
  virtual bool endIteration(Rvector& solution, const Rvector& residual) {
    return endIteration(vectorRefSet(1, solution), constVectorRefSet(1, residual));
  }

  void solution(const std::vector<int>& roots, vectorRefSet parameters, vectorRefSet residual,
                vectorRefSetP parametersP = nullVectorRefSetP<Rvector>) {
    auto working_set_save = m_working_set;
    m_working_set = roots;
    buildSubspace(true);
    solveReducedProblem();
    doInterpolation(parameters, residual, parametersP);
    m_working_set = working_set_save;
  }
  void solution(const std::vector<int>& roots, std::vector<Rvector>& parameters, std::vector<Rvector>& residual,
                vectorSetP& parametersP = nullVectorSetP<Rvector>) {
    return solution(roots, vectorRefSet(parameters.begin(), parameters.end()),
                    vectorRefSet(residual.begin(), residual.end()),
                    vectorRefSetP(parametersP.begin(), parametersP.end()));
  }
  void solution(int root, Rvector& parameters, Rvector& residual, vectorP& parametersP = nullVectorP<Rvector>) {
    return solution(std::vector<int>(1, root), vectorRefSet(1, parameters), vectorRefSet(1, residual),
                    vectorRefSetP(1, parametersP));
  }

  /*!
   * \brief Get the solver's suggestion of which degrees of freedom would be best
   * to add to the P-space.
   * \param solution Current solution
   * \param residual Current residual
   * \param maximumNumber Suggest no more than this number
   * \param threshold Suggest only axes for which the current residual and update
   * indicate an energy improvement in the next iteration of this amount or more.
   * \return
   */
  std::vector<size_t> suggestP(constVectorRefSet solution, constVectorRefSet residual,
                               const size_t maximumNumber = 1000, const double threshold = 0) {
    std::map<size_t, scalar_type> result;
    for (size_t kkk = 0; kkk < solution.size(); kkk++) {
      {
        std::vector<size_t> indices;
        std::vector<scalar_type> values;
        //        std::tie(indices, values) = solution[kkk].get().select(residual[kkk], maximumNumber, threshold);
        auto selection = m_handlers->rr().select_max_dot(maximumNumber, solution[kkk], solution[kkk]);
        indices.reserve(selection.size());
        values.reserve(selection.size());
        for (auto elem : selection) {
          indices.emplace_back(elem.first);
          values.emplace_back(elem.second);
        }
        //             molpro::cout <<"indices.size()="<<indices.size()<<std::endl;
        //             for (auto k=0; k<indices.size(); k++) molpro::cout << "select "<< indices[k] <<" :
        //             "<<values[k]<<std::endl;
        for (size_t i = 0; i < indices.size(); i++)
          if (result.count(indices[i]))
            result[indices[i]] = std::max(result[indices[i]], values[i]);
          else
            result[indices[i]] = values[i];
      }
    }
    // sort and select
    //   for (const auto& kv : result) molpro::cout << "result: " << kv.first << " : " <<kv.second<<std::endl;
    std::multimap<scalar_type, size_t, std::greater<scalar_type>> inverseResult;
    for (const auto& kv : result)
      inverseResult.insert(std::pair<scalar_type, size_t>(kv.second, kv.first));
    //   for (const auto& kv : inverseResult) molpro::cout << "inverseResult: " << kv.first << " : "
    //   <<kv.second<<std::endl;
    std::vector<size_t> indices;
    //   std::vector<T> values;
    size_t k = 0;
    for (auto p = inverseResult.cbegin(); p != inverseResult.cend() && k < maximumNumber; k++) {
      indices.push_back(p->second); // values.push_back(p->first);
      ++p;
    }
    //   for (auto k=0; k<indices.size(); k++) molpro::cout << "suggest P "<< indices[k] <<" : "<<values[k]<<std::endl;
    return indices;
  }
  std::vector<size_t> suggestP(const std::vector<Rvector>& solution, const std::vector<Rvector>& residual,
                               const size_t maximumNumber = 1000, const double threshold = 0) {
    return suggestP(constVectorRefSet(solution.begin(), solution.end()),
                    constVectorRefSet(residual.begin(), residual.end()), maximumNumber, threshold);
  }

  /*!
   * \brief Set convergence threshold
   */
  void setThresholds(double thresh) { m_thresh = thresh; }

  std::vector<scalar_type> eigenvalues() const ///< The calculated eigenvalues of the subspace matrix
  {
    std::vector<scalar_type> result;
    for (size_t root = 0; root < (size_t)m_roots && root < m_eval_xx.size(); root++)
      result.push_back(m_eval_xx[root]);
    return result;
  }

  /*!
   * @brief The roots that are currently being tracked
   * @return
   */
  const std::vector<int>& working_set() const { return m_working_set; }

  std::vector<scalar_type>
  working_set_eigenvalues() const ///< The calculated eigenvalues of the subspace matrix belonging to the working set
  {
    std::vector<scalar_type> result;
    for (const auto& root : m_working_set)
      result.push_back(m_eval_xx[root]);
    return result;
  }

  const std::vector<scalar_type>& errors() const { return m_errors; } //!< Error at last iteration

public:
  int m_verbosity; //!< How much to print. Zero means nothing; One results in a single progress-report line printed each
                   //!< iteration.
  double m_thresh; //!< If residual . residual is less than this, converged.
  unsigned int m_maxIterations; //!< Maximum number of iterations
  unsigned int m_minIterations; //!< Minimum number of iterations
  bool m_linear;                ///< Whether residuals are linear functions of the corresponding expansion vectors.
  bool m_hermitian; ///< Whether residuals can be assumed to be the action of an underlying self-adjoint operator.
  size_t
      m_roots; ///< How many roots to calculate / equations to solve (defaults to size of solution and residual vectors)
  bool m_rspt;

public:
  optionMap m_options; ///< A string of options to be interpreted by solveReducedProblem().
  ///< Possibilities include
  ///< - m_options["convergence"]=="energy" (the default), meaning that m_errors() returns the predicted eigenvalue
  ///< change in the next iteration, ie the scalar product of the step and the residual
  ///< - m_options["convergence"]=="step": m_errors() returns the norm of the step in the solution
  ///< - m_options["convergence"]=="residual": m_errors() returns the norm of the residual vector
protected:
  iterativesolver::Q<Rvector, Qvector, Pvector, scalar_type> m_qspace;
  iterativesolver::P<Pvector> m_pspace;
  std::vector<Qvector> m_last_d;    ///< optimum solution in last iteration
  std::vector<Qvector> m_last_hd;   ///< action vector corresponding to optimum solution in last iteration
  std::vector<Qvector> m_current_r; ///< current working space TODO can probably eliminate using m_last_d
  std::vector<Qvector> m_current_v; ///< action vector corresponding to current working space
  std::vector<std::vector<value_type>> m_s_rr, m_h_rr, m_hh_rr, m_rhs_r;  ///< interactions within R space
  std::map<int, std::vector<value_type>> m_s_qr, m_h_qr, m_h_rq, m_hh_qr; ///< interactions between R and Q spaces
  std::vector<std::vector<value_type>> m_s_pr, m_h_pr, m_h_rp;            ///< interactions between R and P spaces
  mutable std::vector<int> m_working_set; ///< which roots are being tracked in the working set
  std::map<int, int> m_q_solutions;       ///< key of q space vector of a converged solution
  bool m_exclude_r_from_redundancy_test;
  bool m_orthogonalise_Q; //!< whether Q-space vectors constructed by difference should be orthogonal to the working
                          //!< vector, or the pure difference with the previous vector

public:
  virtual void report() const {
    if (m_verbosity > 0) {
      molpro::cout << "iteration " << m_statistics.iterations;
      if (not m_values.empty())
        molpro::cout << ", " << m_value_print_name << " = " << m_values.back();
      if (this->m_roots > 1)
        molpro::cout << ", error[" << std::max_element(m_errors.cbegin(), m_errors.cend()) - m_errors.cbegin()
                     << "] = ";
      else
        molpro::cout << ", error = ";
      molpro::cout << *std::max_element(m_errors.cbegin(), m_errors.cend()) << std::endl;
    }
  }

protected:
  virtual bool solveReducedProblem() = 0;

  void buildSubspace(bool emptyR = false) {
    const size_t nP = m_pspace.size();
    const size_t nQ = m_qspace.size();
    const size_t nR = emptyR ? 0 : m_working_set.size();
    auto& nX = m_n_x;
    nX = nP + nQ + nR;
    const auto oP = 0;
    const auto oQ = oP + nP;
    const auto oR = oQ + nQ;
    //    molpro::cout << "buildSubspace nP=" << nP << ", nQ=" << nQ << ", nR=" << nR << std::endl;
    m_h_xx.resize(nX * nX);
    m_s_xx.resize(nX * nX);
    m_rhs_x.resize(nX * m_rhs.size());
    for (size_t a = 0; a < nQ; a++) {
      for (size_t rhs = 0; rhs < m_rhs.size(); rhs++)
        m_rhs_x[oQ + a + nX * rhs] = m_qspace.rhs(a)[rhs];
      for (size_t b = 0; b < nQ; b++) {
        m_h_xx[oQ + b + nX * (oQ + a)] = m_qspace.action(b, a);
        m_s_xx[oQ + b + nX * (oQ + a)] = m_qspace.metric(b, a);
      }
      const auto& metric_pspace = m_qspace.metric_pspace(a);
      const auto& action_pspace = m_qspace.action_pspace(a);
      for (size_t i = 0; i < nP; i++) {
        m_h_xx[oP + i + nX * (oQ + a)] = m_h_xx[oQ + a + nX * (oP + i)] = action_pspace[i];
        m_s_xx[oP + i + nX * (oQ + a)] = m_s_xx[oQ + a + nX * (oP + i)] = metric_pspace[i];
      }
      for (size_t m = 0; m < nR; m++) {
        m_h_xx[oR + m + nX * (oQ + a)] = m_h_rq[a][m];
        m_h_xx[oQ + a + nX * (oR + m)] = m_h_qr[a][m];
        m_s_xx[oR + m + nX * (oQ + a)] = m_s_qr[a][m];
        m_s_xx[oQ + a + nX * (oR + m)] = m_s_qr[a][m];
      }
    }
    for (size_t i = 0; i < nP; i++) {
      for (size_t m = 0; m < nR; m++) {
        m_h_xx[oR + m + nX * (oP + i)] = m_h_rp[i][m];
        m_h_xx[oP + i + nX * (oR + m)] = m_h_pr[i][m];
        m_s_xx[oR + m + nX * (oP + i)] = m_s_pr[i][m];
        m_s_xx[oP + i + nX * (oR + m)] = m_s_pr[i][m];
      }
      for (size_t j = 0; j < nP; j++) {
        m_h_xx[oP + i + nX * (oP + j)] = m_pspace.action(i, j);
        m_s_xx[oP + i + nX * (oP + j)] = m_pspace.metric(i, j);
      }
    }
    for (size_t n = 0; n < nR; n++) {
      for (size_t rhs = 0; rhs < m_rhs.size(); rhs++)
        m_rhs_x[oR + n + nX * rhs] = m_rhs_r[n][rhs];
      for (size_t m = 0; m < nR; m++) {
        m_h_xx[oR + m + nX * (oR + n)] = m_h_rr[m][n];
        m_s_xx[oR + m + nX * (oR + n)] = m_s_rr[m][n];
      }
    }
    if (m_subspaceMatrixResRes)
      m_s_xx = m_h_xx;
    if (nQ > 0) {
      std::vector<size_t> candidates;
      std::map<int, int> solutions_q; // map from current q space to roots
      for (const auto& x :
           m_q_solutions) // loop over all roots that are associated with a q vector. x.first=root, x.second=qindex
        for (auto k = 0; k < m_qspace.keys().size(); k++)
          if (m_qspace.keys()[k] == x.second)
            solutions_q[k] = x.first; // TODO check logic
      for (auto a = 0; a < nQ; a++)
        if (solutions_q.count(a) == 0)
          candidates.push_back(oQ + a);
      //      molpro::cout << "candidates:";
      //      for (const auto& c : candidates)
      //        molpro::cout << " " << c;
      //      molpro::cout << std::endl;
      for (auto k = 0; k < nX; k++)
        if (std::abs(m_s_xx[k * (nX + 1)] - 1) < 1e-15)
          m_s_xx[k * (nX + 1)] =
              1; // somehow avoid problems that eigen with Intel 18 get the SVD wrong if near-unit matrix
      auto del = iterativesolver::propose_singularity_deletion(m_exclude_r_from_redundancy_test ? nX - nR : nX, nX,
                                                               m_residual_eigen ? m_s_xx.data() : m_h_xx.data(),
                                                               candidates, nQ > m_maxQ ? 1e6 : m_singularity_threshold);
      if (del >= 0) {
        if (m_verbosity > 2)
          molpro::cout << "del=" << del << "; remove Q" << del - oQ << std::endl;
        m_statistics.q_deletions++;
        m_qspace.remove(del - oQ);
        m_errors.assign(m_roots, 1e20);
        for (auto m = 0; m < nR; m++)
          for (auto a = del - oQ; a < nQ - 1; a++) {
            m_h_rq[a][m] = m_h_rq[a + 1][m];
            m_h_qr[a][m] = m_h_qr[a + 1][m];
            m_s_qr[a][m] = m_s_qr[a + 1][m];
            m_hh_qr[a][m] = m_hh_qr[a + 1][m];
          }
        buildSubspace(emptyR);
        return;
      }
    }
    if (m_verbosity > 1)
      molpro::cout << "nP=" << nP << ", nQ=" << nQ << ", nR=" << nR << std::endl;
    if (m_verbosity > 2) {
      iterativesolver::printMatrix(this->m_s_xx, nX, nX, "Subspace overlap");
      iterativesolver::printMatrix(this->m_h_xx, nX, nX, "Subspace matrix");
    }
  }

protected:
  void diagonalizeSubspaceMatrix() {
    iterativesolver::eigenproblem(m_evec_xx, m_eval_xx, m_h_xx, m_s_xx, m_n_x, m_hermitian, m_svdThreshold,
                                  m_verbosity);
  }

  /*!
   * @brief form the combination of P, Q and R vectors to give the interpolated solution and corresponding residual
   *
   * @param solution On exit, the complete current solution (R, P and Q parts)
   * @param residual On exit, the R and Q contribution to the residual. The action of the matrix on the P solution is
   * missing, and has to be evaluated by the caller.
   * @param solutionP On exit, the solution projected to the P space
   * @param actionOnly If true, omit P space contribution and calculate action vector, not full residual
   */
  void doInterpolation(vectorRefSet solution, vectorRefSet residual, vectorRefSetP solutionP,
                       bool actionOnly = false) const {
    for (auto& s : solution)
      m_handlers->rr().fill(0, s);
    for (auto& s : residual)
      m_handlers->rr().fill(0, s);
    const auto nP = m_pspace.size();
    const auto nR = m_current_r.size();
    //    auto nQ = m_qspace.size();
    // FIXME Why is this a reference?
    const auto& nX = this->m_solution_x.size() / this->m_roots;
    const auto nQ =
        nX - nP - nR; // guard against using any vectors added to the Q space since the subspace solution was evaluated
    assert(nQ <= m_qspace.size());
    const size_t oP = 0;
    const auto oQ = oP + nP;
    const auto oR = oQ + nQ;
    assert(m_working_set.size() <= solution.size());
    assert(nP == 0 || solutionP.size() == residual.size());
    for (size_t kkk = 0; kkk < m_working_set.size(); kkk++) {
      auto root = m_working_set[kkk];
      //      molpro::cout << "working set k=" << kkk << " root=" << root << std::endl;
      if (nP > 0)
        solutionP[kkk].get().resize(nP);
      if (not actionOnly) {
        for (size_t i = 0; i < nP; i++) {
          m_handlers->rp().axpy((solutionP[kkk].get()[i] = this->m_solution_x[oP + i + nX * root]), m_pspace[i],
                                solution[kkk]);
        }
      }
      //      molpro::cout << "square norm of solution after P contribution " << solution[kkk]->dot(*solution[kkk]) <<
      //      std::endl;
      for (size_t q = 0; q < nQ; q++) {
        auto l = oQ + q;
        m_handlers->rq().axpy(this->m_solution_x[l + nX * root], m_qspace[q], solution[kkk]);
        m_handlers->rq().axpy(this->m_solution_x[l + nX * root], m_qspace.action(q), residual[kkk]);
      }
      if (true) {
        for (int c = 0; c < nR; c++) {
          auto l = oR + c;
          m_handlers->rq().axpy(this->m_solution_x[l + nX * root], m_current_r[c], solution[kkk]);
          m_handlers->rq().axpy(this->m_solution_x[l + nX * root], m_current_v[c], residual[kkk]);
        }
        if (m_residual_eigen) {
          auto norm = m_handlers->rr().dot(solution[kkk], solution[kkk]);
          if (norm != 0) {
            m_handlers->rr().scal(1 / std::sqrt(norm), solution[kkk]);
            m_handlers->rr().scal(1 / std::sqrt(norm), residual[kkk]);
          }
        }
        // TODO
      }
      if (not actionOnly and (m_residual_eigen || (m_residual_rhs && m_augmented_hessian > 0)))
        m_handlers->rr().axpy(-this->m_eval_xx[root], solution[kkk], residual[kkk]);
      if (not actionOnly and m_residual_rhs)
        m_handlers->rr().axpy(-1, this->m_rhs[root], residual[kkk]);
    }
  }

public:
  const iterativesolver::Statistics& statistics() const { return m_statistics; }

public:
  std::vector<double> m_errors; //!< Error at last iteration
  bool m_subspaceMatrixResRes;  // whether m_subspaceMatrix is Residual.Residual (true) or Solution.Residual (false)
  bool m_residual_eigen;        // whether to subtract eigenvalue*solution when constructing residual
  bool m_residual_rhs;          // whether to subtract rhs when constructing residual
  // whether to use RSPT to construct solution instead of diagonalisation
  vectorSet m_rhs;
  size_t m_lastVectorIndex;
  std::vector<scalar_type> m_updateShift;
  size_t m_n_x;                         //!< size of full subspace
  std::vector<value_type> m_h_xx;       //!< full subspace
  std::vector<value_type> m_s_xx;       //!< full subspace
  std::vector<value_type> m_rhs_x;      //!< full subspace
  std::vector<value_type> m_evec_xx;    //!< full subspace
  std::vector<value_type> m_eval_xx;    //!< full subspace
  std::vector<value_type> m_values;     //!< function values
  std::vector<value_type> m_solution_x; //!< solution in x space
public:
  size_t m_dimension;             //!< not used in the class, but a place for clients (eg C interface) to store a number
                                  //!< representing the size of the underlying vector space.
  std::string m_value_print_name; //< the title report() will give to the function value
protected:
  double m_singularity_threshold;
  double m_augmented_hessian; //!< The scale factor for augmented hessian solution of linear inhomogeneous systems.
                              //!< Special values:
                              //!< - 0: unmodified linear equations
                              //!< - 1: standard augmented hessian

  bool m_nullify_solution_before_update;
  iterativesolver::Statistics m_statistics;

public:
  double m_svdThreshold; ///< Threshold for singular-value truncation in linear equation solver.
  size_t m_maxQ;         //!< maximum size of Q space
protected:
};

/*! @example LinearEigensystemExample.cpp */
/*!
 * \brief A class that finds the lowest eigensolutions of a matrix using Davidson's method, i.e. preconditioned Lanczos
 *
 * Example of simplest use with a simple in-memory container for eigenvectors: @include LinearEigensystemExample.cpp
 *
 * Example using a P-space and offline distributed storage provided by the PagedVector class: @include
 * LinearEigensystemExample-paged.cpp
 *
 * @tparam Rvector The class encapsulating solution and residual vectors
 * @tparam Qvector Used internally as a class for storing vectors on backing store
 * @tparam Pvector Specify a P-space vector as a sparse combination of parameters.
 */
template <class Rvector = std::vector<double>, class Qvector = Rvector,
          class Pvector = std::map<size_t, typename Rvector::value_type>>
class LinearEigensystem : public IterativeSolver<Rvector, Qvector, Pvector> {
public:
  using typename IterativeSolver<Rvector, Qvector, Pvector>::scalar_type;
  using typename IterativeSolver<Rvector, Qvector, Pvector>::value_type;
  using IterativeSolver<Rvector, Qvector, Pvector>::m_verbosity;

  /*!
   * \brief LinearEigensystem
   */
  explicit LinearEigensystem(const std::shared_ptr<iterativesolver::ArrayHandlers<Rvector, Qvector, Pvector>>& handlers)
      : IterativeSolver<Rvector, Qvector, Pvector>(handlers) {
    this->m_residual_rhs = false;
    this->m_residual_eigen = true;
    this->m_linear = true;
  }

private:
  bool solveReducedProblem() override {
    if (this->m_rspt) {
      throw std::logic_error("RSPT not yet implemented");
    } else {
      this->diagonalizeSubspaceMatrix();
      //                   << std::endl;
      this->m_solution_x.resize(this->m_n_x * std::min(int(this->m_roots), int(this->m_n_x)));
      std::copy_n(this->m_evec_xx.begin(), this->m_solution_x.size(), this->m_solution_x.begin());
    }

    this->m_updateShift.resize(this->m_roots);
    for (size_t root = 0; root < (size_t)this->m_roots; root++)
      this->m_updateShift[root] =
          -(1 + std::numeric_limits<scalar_type>::epsilon()) * (root < this->m_n_x ? this->m_eval_xx[root] : 0);
    return true;
  }

public:
  void report() const override {
    std::vector<value_type> ev = this->eigenvalues();
    if (this->m_verbosity > 0) {
      molpro::cout << "iteration " << this->m_statistics.iterations << "[" << this->m_working_set.size() << "]";
      if (this->m_pspace.size() > 0)
        molpro::cout << ", P=" << this->m_pspace.size();
      if (this->m_roots > 1)
        molpro::cout << ", error["
                     << std::max_element(this->m_errors.cbegin(), this->m_errors.cend()) - this->m_errors.cbegin()
                     << "] = ";
      else
        molpro::cout << ", error = ";
      molpro::cout << *std::max_element(this->m_errors.cbegin(), this->m_errors.cend()) << ", eigenvalues: ";
      for (const auto e : ev)
        molpro::cout << " " << e;
      molpro::cout << std::endl;
    }
  }
};

/** @example LinearEquationsExample.cpp */
/*!
 * \brief A class that finds the solutions of linear equation systems using a generalisation of Davidson's method, i.e.
 * preconditioned Lanczos
 *
 * Example of simplest use: @include LinearEquationsExample.cpp
 * @tparam Rvector The class encapsulating solution and residual vectors
 * @tparam Qvector Used internally as a class for storing vectors on backing store
 * @tparam Pvector Specify a P-space vector as a sparse combination of parameters.
 */
template <class Rvector = std::vector<double>, class Qvector = Rvector,
          class Pvector = std::map<size_t, typename Rvector::value_type>>
class LinearEquations : public IterativeSolver<Rvector, Qvector, Pvector> {
protected:
  using typename IterativeSolver<Rvector, Qvector, Pvector>::scalar_type;
  using typename IterativeSolver<Rvector, Qvector, Pvector>::value_type;

public:
  using vectorSet = typename std::vector<Rvector>;                                       ///< Container of vectors
  using vectorRefSet = typename std::vector<std::reference_wrapper<Rvector>>;            ///< Container of vectors
  using constVectorRefSet = typename std::vector<std::reference_wrapper<const Rvector>>; ///< Container of vectors
  using IterativeSolver<Rvector, Qvector, Pvector>::m_verbosity;

  /*!
   * \brief Constructor
   * \param rhs right-hand-side vectors. More can be added subsequently using addEquations(), provided iterations have
   * not yet started. \param augmented_hessian If zero, solve the inhomogeneous equations unmodified. If 1, solve
   * instead the augmented hessian problem. Other values scale the augmented hessian damping.
   * \param handlers group of array handlers for coordinating array operations
   */
  explicit LinearEquations(constVectorRefSet rhs,
                           const std::shared_ptr<iterativesolver::ArrayHandlers<Rvector, Qvector, Pvector>>& handlers,
                           double augmented_hessian = 0)
      : IterativeSolver<Rvector, Qvector, Pvector>(handlers) {
    this->m_linear = true;
    this->m_residual_rhs = true;
    this->m_augmented_hessian = augmented_hessian;
    addEquations(rhs);
  }

  explicit LinearEquations(const vectorSet& rhs,
                           const std::shared_ptr<iterativesolver::ArrayHandlers<Rvector, Qvector, Pvector>>& handlers,
                           double augmented_hessian = 0)
      : LinearEquations(constVectorRefSet(rhs.begin(), rhs.end()), handlers, augmented_hessian) {}

  explicit LinearEquations(const Rvector& rhs,
                           const std::shared_ptr<iterativesolver::ArrayHandlers<Rvector, Qvector, Pvector>>& handlers,
                           double augmented_hessian = 0)
      : LinearEquations(constVectorRefSet(1, rhs), handlers, augmented_hessian) {}

  /*!
   * \brief add one or more equations to the set to be solved, by specifying their right-hand-side vector
   * \param rhs right-hand-side vectors to be added
   */
  void addEquations(constVectorRefSet rhs) {
    //   for (const auto &v : rhs) this->m_rhs.push_back(v);
    this->m_rhs.clear();
    this->m_rhs.reserve(rhs.size());
    for (const auto& v : rhs)
      // FIXME Is this meant to be a copy?
      this->m_rhs.emplace_back(this->m_handlers->rr().copy(v)); // TODO template-ise these options
    //   molpro::cout << "addEquations makes m_rhs.back()="<<this->m_rhs.back()<<std::endl;
  }
  void addEquations(const std::vector<Rvector>& rhs) { addEquations(vectorSet(rhs.begin(), rhs.end())); }
  void addEquations(const Rvector& rhs) { addEquations(vectorSet(1, rhs)); }

protected:
  bool solveReducedProblem() override {
    iterativesolver::solve_LinearEquations(this->m_solution_x, this->m_eval_xx, this->m_h_xx, this->m_s_xx,
                                           this->m_rhs_x, this->m_n_x, this->m_roots, this->m_augmented_hessian,
                                           this->m_svdThreshold, this->m_verbosity);
    return true;
  }
};

/** @example OptimizeExample.cpp */
/*!
 * \brief A class that optimises a function using a Quasi-Newton or other method
 *
 * Example of simplest use: @include OptimizeExample.cpp
 * @tparam Rvector The class encapsulating solution and residual vectors
 * @tparam Qvector Used internally as a class for storing vectors on backing store
 */
template <class Rvector = std::vector<double>, class Qvector = Rvector>
class Optimize : public IterativeSolver<Rvector, Qvector> {
public:
  using typename IterativeSolver<Rvector, Qvector>::scalar_type;
  using typename IterativeSolver<Rvector, Qvector>::value_type;
  using vectorSet = typename std::vector<Rvector>;                                       ///< Container of vectors
  using vectorRefSet = typename std::vector<std::reference_wrapper<Rvector>>;            ///< Container of vectors
  using constVectorRefSet = typename std::vector<std::reference_wrapper<const Rvector>>; ///< Container of vectors
  using IterativeSolver<Rvector, Qvector>::m_verbosity;
  using IterativeSolver<Rvector, Qvector>::m_values;

  /*!
   * \brief Constructor
   * \param algorithm Allowed values: "L-BFGS","null"
   * \param minimize If false, a maximum, not minimum, will be sought
   * \param handlers group of array handlers for coordinating array operations
   */
  explicit Optimize(
      const std::shared_ptr<iterativesolver::ArrayHandlers<Rvector, Qvector, std::map<size_t, double>>>& handlers,
      std::string algorithm = "L-BFGS", bool minimize = true)
      : IterativeSolver<Rvector, Qvector>(handlers), m_algorithm(std::move(algorithm)), m_minimize(minimize),
        m_strong_Wolfe(true), m_Wolfe_1(0.0001), m_Wolfe_2(0.9), // recommended values Nocedal and Wright p142
        m_linesearch_tolerance(0.2), m_linesearch_grow_factor(3), m_linesearch_steplength(0) {
    this->m_linear = false;
    this->m_residual_rhs = false;
    this->m_residual_eigen = false;
    this->m_roots = 1;
    this->m_subspaceMatrixResRes = false;
    this->m_singularity_threshold = 0;
    this->m_orthogonalise_Q = false;
    this->m_exclude_r_from_redundancy_test = true;
    this->m_hermitian = false;
    this->m_qspace.hermitian(false);
  }

protected:
  std::string m_algorithm; ///< which variant of Quasi-Newton or other methods
  bool m_minimize;         ///< whether to minimize or maximize
  using IterativeSolver<Rvector, Qvector>::m_handlers;

public:
  bool m_strong_Wolfe;             /// Whether to use strong or weak Wolfe conditions
  double m_Wolfe_1;                ///< Acceptance parameter for function value
  double m_Wolfe_2;                ///< Acceptance parameter for function gradient
  double m_linesearch_tolerance;   ///< If the predicted line search is within tolerance of 1, don't bother taking it
  double m_linesearch_grow_factor; ///< If the predicted line search step is extrapolation, limit the step to this
                                   ///< factor times the current step
protected:
  double m_linesearch_steplength; ///< the current line search step. Zero means continue with QN
  //  double m_linesearch_quasinewton_steplength; ///< what fraction of the Quasi-Newton step is the current line
  //  search step
  std::unique_ptr<Qvector> m_best_r, m_best_v;
  double m_best_f;

  bool interpolatedMinimum(value_type& x, double& f, value_type x0, value_type x1, double f0, double f1, value_type g0,
                           value_type g1) {
    if (std::abs(2 * f1 - g1 - 2 * f0 - g0) < 1e-10) { // cubic coefficient is zero
      auto c2 = (g1 - g0) / 2;
      if (c2 < 0)
        return false;
      x = x0 + (-0.5 * g0 / c2) * (x1 - x0);
      f = f0 + g0 * x + c2 * x * x;
      return true;
    }
    auto discriminant = (std::pow(3 * f0 - 3 * f1 + g0, 2) + (6 * f0 - 6 * f1 + g0) * g1 + std::pow(g1, 2));
    //    molpro::cout << "discriminant " << discriminant << std::endl;
    if (discriminant < 0)
      return false; // cubic has no turning points

    auto alpham = (2 * f0 - 2 * f1 + g0 + g1 == 0)
                      ? (g0 / (2 * f1 - 2 * f0 - 2 * g1))
                      : (3 * f0 - 3 * f1 + 2 * g0 + g1 - std::sqrt(discriminant)) / (3 * (2 * f0 - 2 * f1 + g0 + g1));
    auto alphap = (2 * f0 - 2 * f1 + g0 + g1 == 0)
                      ? (g0 / (2 * f1 - 2 * f0 - 2 * g1))
                      : (3 * f0 - 3 * f1 + 2 * g0 + g1 + std::sqrt(discriminant)) / (3 * (2 * f0 - 2 * f1 + g0 + g1));
    auto fm = f0 + alpham * (g0 + alpham * (-3 * f0 + 3 * f1 - 2 * g0 - g1 + alpham * (2 * f0 - 2 * f1 + g0 + g1)));
    auto fp = f0 + alphap * (g0 + alphap * (-3 * f0 + 3 * f1 - 2 * g0 - g1 + alphap * (2 * f0 - 2 * f1 + g0 + g1)));
    f = std::min(fm, fp);
    x = x0 + (fm < fp ? alpham : alphap) * (x1 - x0);
    return true;
  }

  bool solveReducedProblem() override {
    auto n = this->m_qspace.size();
    assert(this->m_n_x == n + 1);
    this->m_solution_x.assign(n + 1, 0);
    this->m_solution_x.back() = 1;
    //    molpro::cout << "Optimize::solveReduced Problem n=" << n << std::endl;
    if (n > 0) {

      // first consider whether this point can be taken as the next iteration point, or whether further line-searching
      // is needed
      //      auto step = std::sqrt(this->m_subspaceOverlap(n - 1, n - 1));
      double step = 1 / this->m_qspace.scale_factor(this->m_qspace.size() - 1);
      auto f0 = m_best_f;
      auto f1 = this->m_values.back();
      auto g1 = step * this->m_h_qr[n - 1][0];
      auto g0 = step * m_handlers->qq().dot((*m_best_v), this->m_qspace[this->m_qspace.size() - 1]);
      bool Wolfe_1 = f1 <= f0 + m_Wolfe_1 * g0;
      bool Wolfe_2 = m_strong_Wolfe ? g1 >= m_Wolfe_2 * g0 : std::abs(g1) <= m_Wolfe_2 * std::abs(g0);
      if (this->m_verbosity > 1) {
        molpro::cout << "step=" << step << std::endl;
        molpro::cout << "f0=" << f0 << std::endl;
        molpro::cout << "f1=" << f1 << std::endl;
        molpro::cout << " m_Wolfe_1 =" << m_Wolfe_1 << std::endl;
        molpro::cout << " m_Wolfe_1 * g0=" << m_Wolfe_1 * g0 << std::endl;
        molpro::cout << "f0 + m_Wolfe_1 * g0=" << f0 + m_Wolfe_1 * g0 << std::endl;
        molpro::cout << "g0=" << g0 << std::endl;
        molpro::cout << "g1=" << g1 << std::endl;
        molpro::cout << "Wolfe conditions: " << Wolfe_1 << Wolfe_2 << std::endl;
      }
      if (g1 < this->m_thresh or (Wolfe_1 && Wolfe_2))
        goto accept;
      double finterp;
      //      molpro::cout << "before interpolatedMinimum" << std::endl;
      value_type alpha;
      auto interpolated = interpolatedMinimum(alpha, finterp, 0, 1, f0, f1, g0, g1);
      //      molpro::cout << "interpolated: " << interpolated << ", alpha " <<
      //      alpha
      //                   << ", finterp " << finterp << std::endl;
      if (interpolated and ((g0 > 0 and g1 > 0 and alpha > 0) or
                            (g0 < 0 and g1 < 0 and alpha < 1))) // not bracketed, interpolant goes the wrong way
        interpolated = false;
      if (not interpolated or alpha > m_linesearch_grow_factor) {
        if (this->m_verbosity > 1) {
          if (interpolated)
            molpro::cout << "reject interpolated minimum value " << finterp << " at alpha=" << alpha << std::endl;
          else
            molpro::cout << "cubic interpolation did not find a valid minimum" << std::endl;
          molpro::cout << "taking instead step=" << (g0 > 0 ? 1 : -1) * m_linesearch_grow_factor << std::endl;
        }
        alpha = (g0 > 0 ? 1 : -1) * m_linesearch_grow_factor; // expand the search range
      } else if (std::abs(alpha - 1) < m_linesearch_tolerance) {
        if (this->m_verbosity > 1)
          molpro::cout << "Don't bother with linesearch " << alpha << std::endl;
        goto accept; // if we are within spitting distance already, don't bother to make a line step
      } else {
        if (this->m_verbosity > 1)
          molpro::cout << "cubic linesearch interpolant has minimum " << finterp << " at " << alpha << "(absolute step "
                       << (alpha - 1) * step << ")" << std::endl;
      }
      // when we arrive here, we need to do a new line-search step
      //      molpro::cout << "we need to do a new line-search step " << alpha << std::endl;
      if (f1 <= f0) {
        m_linesearch_steplength = (alpha - 1) * step;
        // FIXME make_unique using a copy. Can this be a weak ptr or a reference?
        this->m_statistics.best_r_creations++;
        m_best_r.reset(new Qvector(this->m_current_r.front()));
        m_best_v.reset(new Qvector(this->m_current_v.front()));
        m_best_f = this->m_values.back();
        //        molpro::cout << "setting best to current, with f=" << m_best_f << std::endl;
      } else
        m_linesearch_steplength = alpha * step;

      this->m_nullify_solution_before_update = false;
      return false;
    }
  accept:
    //    molpro::cout << "accept reached" << std::endl;
    m_linesearch_steplength = 0;
    this->m_nullify_solution_before_update = true;
    if (this->m_algorithm == "L-BFGS") {
      for (int aint = this->m_qspace.size() - 1; aint >= 0; aint--) {
        size_t a = std::abs(aint);
        //        molpro::cout << "iterate q_" << a << std::endl;
        this->m_solution_x[a] = -this->m_h_qr[a][0];
        for (auto b = a + 1; b < this->m_qspace.size(); b++)
          this->m_solution_x[a] -= this->m_solution_x[b] * this->m_qspace.action(a, b);
        this->m_solution_x[a] /= this->m_qspace.action(a, a);
      }
    }
    this->m_statistics.best_r_creations++;
    m_best_r.reset(new Qvector(this->m_current_r.front()));
    m_best_v.reset(new Qvector(this->m_current_v.front()));
    m_best_f = this->m_values.back();

    return true;
  }

public:
  virtual bool endIteration(vectorRefSet solution, constVectorRefSet residual) override {
    if (this->m_q_solutions.count(0) == 0) {
      if (m_linesearch_steplength != 0) { // line search
        //        molpro::cout << "*enter endIteration m_linesearch_steplength=" << m_linesearch_steplength <<
        //        std::endl;
        //              molpro::cout << "solution " << solution.front().get() << std::endl;
        // FIXME is this meant to be a copy?
        solution.front().get() = *m_best_r;
        m_handlers->rq().axpy(m_linesearch_steplength, this->m_qspace[this->m_qspace.size() - 1], solution.front());
        this->m_values.pop_back();
        this->m_qspace.remove(this->m_qspace.size() - 1);
      } else { // quasi-Newton
        if (m_algorithm == "L-BFGS" and this->m_solution_x.size() > 0) {
          //          molpro::cout << "L-BFGS stage 2" << std::endl;
          //          molpro::cout << "before subtracting rk solution length="
          //                       << std::sqrt(solution.back().get().dot(solution.back().get())) << std::endl;
          //          solution.back().get().axpy(-1, this->m_last_d.back());
          //          molpro::cout << "after subtracting rk solution length="
          //                       << std::sqrt(solution.back().get().dot(solution.back().get())) << std::endl;
          for (size_t a = 0; a < this->m_qspace.size(); a++) {
            //            molpro::cout << "iterate q_" << a << std::endl;
            auto factor = this->m_solution_x[a] -
                          m_handlers->qr().dot(this->m_qspace.action(a), solution.back()) / this->m_qspace.action(a, a);
            m_handlers->rq().axpy(factor, this->m_qspace[a], solution.back());
            //            molpro::cout << "Q factor " << factor << std::endl;
          }
          //          molpro::cout << "after Q loop solution length=" <<
          //          std::sqrt(solution.back().get().dot(solution.back().get()))
          //                       << std::endl;
          m_handlers->rq().axpy(1, *(this->m_best_r), solution.back());
          //          molpro::cout << "after adding rk solution length="
          //                       << std::sqrt(solution.back().get().dot(solution.back().get())) << std::endl;
        }
      }
      //    molpro::cout << "*exit endIteration m_linesearch_steplength=" << m_linesearch_steplength << std::endl;
    }
    return IterativeSolver<Rvector, Qvector>::endIteration(solution, residual);
  }

  virtual bool endIteration(std::vector<Rvector>& solution, const std::vector<Rvector>& residual) override {
    return endIteration(vectorRefSet(solution.begin(), solution.end()),
                        constVectorRefSet(residual.begin(), residual.end()));
  }
  virtual bool endIteration(Rvector& solution, const Rvector& residual) override {
    return endIteration(vectorRefSet(1, solution), constVectorRefSet(1, residual));
  }

public:
  virtual void report() const override {
    if (this->m_verbosity > 0) {
      molpro::cout << "iteration " << this->m_statistics.iterations;
      if (m_linesearch_steplength != 0)
        molpro::cout << ", line search step = " << m_linesearch_steplength;
      if (not this->m_values.empty())
        molpro::cout << ", " << this->m_value_print_name << " = " << this->m_values.back();
      molpro::cout << ", error = " << this->m_errors.front() << std::endl;
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
 * @tparam Rvector The class encapsulating solution and residual vectors
 * @tparam Qvector Used internally as a class for storing vectors on backing store
 */
template <class Rvector = std::vector<double>, class Qvector = Rvector>
class DIIS : public IterativeSolver<Rvector, Qvector> {

public:
  using typename IterativeSolver<Rvector, Qvector>::scalar_type;
  using typename IterativeSolver<Rvector, Qvector>::value_type;
  using IterativeSolver<Rvector, Qvector>::m_verbosity;

  /*!
   * \brief DIIS
   */
  DIIS(const std::shared_ptr<iterativesolver::ArrayHandlers<Rvector, Qvector, std::map<size_t, double>>>& handlers)
      : IterativeSolver<Rvector, Qvector>(handlers) {
    this->m_residual_rhs = false;
    this->m_residual_eigen = false;
    this->m_roots = 1;
    this->m_exclude_r_from_redundancy_test = true;
    this->m_singularity_threshold =
        this->m_svdThreshold; // It does not matter if the submatrix goes a bit singular in DIIS
    this->m_orthogonalise_Q = false;
  }

protected:
  bool solveReducedProblem() override {
    this->m_updateShift.clear();
    this->m_updateShift.push_back(-(1 + std::numeric_limits<double>::epsilon()) *
                                  this->m_h_xx[0]); // TODO check that this is what is really wanted

    if (this->m_roots > 1)
      throw std::logic_error("DIIS does not handle multiple solutions");

    iterativesolver::solve_DIIS(this->m_solution_x, this->m_h_xx, this->m_n_x, this->m_svdThreshold, this->m_verbosity);
    return true;
  }
};

} // namespace linalg
} // namespace molpro

#endif // ITERATIVESOLVER_H
