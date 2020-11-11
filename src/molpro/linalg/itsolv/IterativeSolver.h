#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_ITERATIVESOLVER_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_ITERATIVESOLVER_H
#include <molpro/linalg/array/Span.h>
#include <molpro/linalg/itsolv/ArrayHandlers.h>
#include <molpro/linalg/itsolv/Statistics.h>
#include <molpro/linalg/itsolv/wrap.h>
#include <ostream>
#include <vector>

namespace molpro::linalg::itsolv {

/*!
 * @brief Base class defining the interface common to all iterative solvers
 *
 * @tparam R container for "working-set" vectors. These are typically implemented in memory, and are created by the
 * client program. R vectors are never created inside IterativeSolver.
 * @tparam Q container for other vectors. These are typically implemented on backing store and/or distributed across
 * processors.  IterativeSolver constructs a number of instances of Q containers to store history.
 * @tparam P a class that specifies the definition of a single P-space vector, which is a strictly sparse vector in the
 * underlying space.
 */
template <class R, class Q, class P>
class IterativeSolver {
public:
  using value_type = typename R::value_type;                          ///< The underlying type of elements of vectors
  using scalar_type = typename array::ArrayHandler<R, Q>::value_type; ///< The type of scalar products of vectors
  using VectorP = std::vector<value_type>; //!< type for vectors projected on to P space, each element is a coefficient
                                           //!< for the corresponding P space parameter
  //! Function type for applying matrix to the P space vectors and accumulating result in a residual
  using fapply_on_p_type = std::function<void(const std::vector<VectorP>&, const CVecRef<P>&, const VecRef<R>&)>;

  virtual ~IterativeSolver() = default;
  IterativeSolver() = default;
  IterativeSolver(const IterativeSolver<R, Q, P>&) = delete;
  IterativeSolver<R, Q, P>& operator=(const IterativeSolver<R, Q, P>&) = delete;
  IterativeSolver(IterativeSolver<R, Q, P>&&) noexcept = default;
  IterativeSolver<R, Q, P>& operator=(IterativeSolver<R, Q, P>&&) noexcept = default;

  virtual size_t add_vector(const VecRef<R>& parameters, const VecRef<R>& action, fapply_on_p_type& apply_p) = 0;
  virtual size_t add_vector(const VecRef<R>& parameters, const VecRef<R>& action, std::vector<VectorP>& pparams) = 0;
  virtual size_t add_vector(const VecRef<R>& parameters, const VecRef<R>& action) = 0;

  // FIXME this should be removed in favour of VecRef interface
  virtual size_t add_vector(std::vector<R>& parameters, std::vector<R>& action, fapply_on_p_type& apply_p) = 0;
  virtual size_t add_vector(std::vector<R>& parameters, std::vector<R>& action, std::vector<VectorP>& pparams) = 0;
  virtual size_t add_vector(std::vector<R>& parameters, std::vector<R>& action) = 0;
  virtual size_t add_vector(R& parameters, R& action) = 0;

  /*!
   * \brief Add P-space vectors to the expansion set for linear methods.
   * \param Pparams the vectors to add. Each Pvector specifies a sparse vector in the underlying space
   * \param pp_action_matrix Matrix projected onto the existing+new, new P space. It should be provided as a
   * 1-dimensional array, with the existing+new index running fastest.
   * \param parameters Used as scratch working space
   * \param action  On exit, the  residual of the interpolated solution.
   * The contribution from the new, and any existing, P parameters is missing, and should be added in subsequently.
   * \param parametersP On exit, the interpolated solution projected onto the P space.
   * \return The number of vectors contained in parameters, action, parametersP
   */
  virtual size_t add_p(const CVecRef<P>& pparams, const array::Span<value_type>& pp_action_matrix,
                       const VecRef<R>& parameters, const VecRef<R>& action, std::vector<VectorP>& parametersP) = 0;
  //! Construct solution and residual for a given set of roots
  virtual void solution(const std::vector<int>& roots, const VecRef<R>& parameters, const VecRef<R>& residual) = 0;
  //! Constructs parameters of selected roots
  virtual void solution_params(const std::vector<int>& roots, const VecRef<R>& parameters) = 0;
  //! Behaviour depends on the solver
  virtual size_t end_iteration(const VecRef<R>& parameters, const VecRef<R>& residual) = 0;
  virtual std::vector<size_t> suggest_p(const CVecRef<R>& solution, const CVecRef<R>& residual, size_t maximumNumber,
                                        double threshold) = 0;

  virtual void solution(const std::vector<int>& roots, std::vector<R>& parameters, std::vector<R>& residual) = 0;
  virtual void solution_params(const std::vector<int>& roots, std::vector<R>& parameters) = 0;
  virtual size_t end_iteration(std::vector<R>& parameters, std::vector<R>& action) = 0;

  /*!
   * @brief Working set of roots that are not yet converged
   */
  virtual const std::vector<int>& working_set() const = 0;
  //! Total number of roots we are solving for, including the ones that are already converged
  virtual size_t n_roots() const = 0;
  virtual void set_n_roots(size_t nroots) = 0;
  virtual const std::vector<scalar_type>& errors() const = 0;
  virtual const Statistics& statistics() const = 0;
  //! Writes a report to cout output stream
  virtual void report(std::ostream& cout) const = 0;
  //! Writes a report to std::cout
  virtual void report() const = 0;

  //! Sets the convergence threshold
  virtual void set_convergence_threshold(double thresh) = 0;
};

/*!
 * @brief Interface for a specific iterative solver, it can add special member functions or variables.
 */
template <class R, class Q, class P>
class LinearEigensystem : public IterativeSolver<R, Q, P> {
public:
  using typename IterativeSolver<R, Q, P>::scalar_type;
  virtual std::vector<scalar_type> eigenvalues() const = 0; //!< The calculated eigenvalues of the subspace matrix
};

template <class R, class Q, class P>
class LinearEquations : public IterativeSolver<R, Q, P> {
public:
  using typename IterativeSolver<R, Q, P>::scalar_type;
  //! eigenvalues of augmented Hessian method, if it was used
  virtual std::vector<scalar_type> eigenvalues() const = 0;
  void addEquations(const std::vector<R>& rhs) = 0;
  std::vector<Q>& rhs() = 0;

protected:
  std::vector<Q> m_rhs;
};

template <class R, class Q, class P>
class Optimize : public IterativeSolver<R, Q, P> {
public:
  virtual bool end_iteration(std::vector<R>& solution, std::vector<R>& residual) = 0;
};

template <class R, class Q, class P>
class DIIS : public IterativeSolver<R, Q, P> {};

} // namespace molpro::linalg::itsolv

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_ITERATIVESOLVER_H
