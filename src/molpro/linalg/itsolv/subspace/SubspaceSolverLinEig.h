#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SUBSPACESOLVERLINEIG_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SUBSPACESOLVERLINEIG_H
#include <molpro/linalg/itsolv/subspace/ISubspaceSolver.h>
#include <molpro/linalg/itsolv/subspace/IXSpace.h>
#include <molpro/linalg/itsolv/subspace/Matrix.h>
#include <molpro/linalg/itsolv/Logger.h>
#include <molpro/linalg/itsolv/helper.h>

#include <map>
#include <memory>
#include <optional>

namespace molpro::linalg::itsolv::subspace {

/*!
 * @brief Solves subspace problem for linear eigenvalues and system of linear equations
 */
template <class RT, class QT, class PT>
class SubspaceSolverLinEig : public ISubspaceSolver<RT, QT, PT> {
public:
  using value_type = typename ISubspaceSolver<RT, QT, PT>::value_type;
  using value_type_abs = typename ISubspaceSolver<RT, QT, PT>::value_type_abs;
  using R = typename ISubspaceSolver<RT, QT, PT>::R;
  using Q = typename ISubspaceSolver<RT, QT, PT>::Q;
  using P = typename ISubspaceSolver<RT, QT, PT>::P;

  explicit SubspaceSolverLinEig(std::shared_ptr<Logger> logger) : m_logger(std::move(logger)) {}

  void solve(IXSpace<R, Q, P>& xspace, const size_t nroots_max) override {
    m_logger->trace("SubspaceSolverLinEig::solve");
    if (xspace.data.find(EqnData::rhs) == xspace.data.end() || xspace.data[EqnData::rhs].empty()) {
      solve_eigenvalue(xspace, nroots_max);
    } else {
      solve_linear_equations(xspace);
    }
  }

protected:
  int convert_verbosity(log::Verbosity verbosity) {
    switch (verbosity) {
      case log::Verbosity::Trace:
        return 3;
      case log::Verbosity::Debug:
        return 2;
      case log::Verbosity::Info:
        return 1;
      case log::Verbosity::None:
        return 0;
    }

    return 0;
  }

  void solve_eigenvalue(IXSpace<R, Q, P>& xspace, const size_t nroots_max) {
    m_logger->trace("SubspaceSolverLinEig::solve_eigenvalue");
    auto h = xspace.data[EqnData::H];
    auto s = xspace.data[EqnData::S];
    m_logger->data_dump("S = ", s);
    m_logger->data_dump<15>("H = ", h);
    auto dim = h.rows();
    auto evec = std::vector<value_type>{};
    int verbosity = convert_verbosity(m_logger->verbosity());
    itsolv::eigenproblem(evec, m_eigenvalues, h.data(), s.data(), dim, m_hermitian, m_svd_solver_threshold, verbosity,
                         &m_imag_eigval_comps, nroots_max);
    m_complex_partners.clear();
    for (size_t k = 0; k + 1 < m_imag_eigval_comps.size(); k += 2) {
      const auto first = m_imag_eigval_comps[k].first;
      const auto second = m_imag_eigval_comps[k + 1].first;
      m_complex_partners[first] = second;
      m_complex_partners[second] = first;
    }
    size_t n_solutions = 0;
    if (dim)
      n_solutions = evec.size() / dim;
    auto full_matrix = Matrix<value_type>{std::move(evec), {n_solutions, dim}};
    auto nroots = std::min(nroots_max, n_solutions);
    m_retained_solutions.clear();
    if (nroots > 0) {
      // Without its partner the residual of a pair member never vanishes. The partner's coefficients are kept, but not
      // the partner itself, so that the number of solutions never exceeds the number of roots requested.
      const auto partner = complex_partner(nroots - 1);
      if (partner && *partner == nroots && *partner < n_solutions) {
        m_retained_solutions.resize({nroots + 1, dim});
        m_retained_solutions.slice() = full_matrix.slice({0, 0}, {nroots + 1, dim});
      }
    }
    m_eigenvalues.resize(nroots);
    auto [first, last] = std::ranges::remove_if(
        m_imag_eigval_comps, [nroots](std::size_t idx) { return idx >= nroots; },
        &decltype(m_imag_eigval_comps)::value_type::first);
    if (first != last && (std::ranges::distance(first, last) % 2) != 0) {
      m_logger->info("Complex eigenvalue pair split up due to truncation to requested number of roots");
    }
    m_imag_eigval_comps.erase(first, last);
    std::erase_if(m_complex_partners, [nroots](const auto& entry) { return entry.first >= nroots; });
    m_solutions.resize({nroots, dim});
    m_solutions.slice() = full_matrix.slice({0, 0}, {nroots, dim});
    m_errors.assign(size(), std::numeric_limits<value_type_abs>::max());
    m_logger->data_dump<10>("eigenvalues = ", m_eigenvalues);
    if (!m_imag_eigval_comps.empty()) {
      m_logger->info("The following eigenvalues turned out to be complex-valued: ",
                     m_imag_eigval_comps |
                         std::ranges::views::transform([](const auto& pair) { return pair.first + 1; }));
      m_logger->data_dump("imaginary parts of eigenvalues = ",
                          m_imag_eigval_comps | std::ranges::views::transform([](auto pair) {
                            pair.first += 1;
                            return pair;
                          }));
    }
    m_logger->data_dump("eigenvectors = ", m_solutions);
  }

  void solve_linear_equations(IXSpace<R, Q, P>& xspace) {
    m_logger->trace("SubspaceSolverLinEig::solve_linear_equations");
    auto h = xspace.data[EqnData::H];
    auto s = xspace.data[EqnData::S];
    auto rhs = xspace.data[EqnData::rhs];
    m_logger->data_dump<15>("S = ", s);
    m_logger->data_dump<15>("H = ", h);
    m_logger->data_dump<15>("rhs = ", rhs);
    const auto dim = h.rows();
    const auto n_solutions = rhs.cols();
    auto solution = std::vector<value_type>{};
    m_eigenvalues.assign(n_solutions, 0);
    m_imag_eigval_comps.clear();
    m_complex_partners.clear();
    m_retained_solutions.clear();
    int verbosity = convert_verbosity(m_logger->verbosity());
    itsolv::solve_LinearEquations(solution, m_eigenvalues, h.data(), s.data(), rhs.data(), dim, n_solutions,
                                  m_augmented_hessian, m_svd_solver_threshold, verbosity);
    m_solutions = Matrix<value_type>{std::move(solution), {n_solutions, dim}};
    m_errors.assign(size(), std::numeric_limits<value_type_abs>::max());
    m_logger->data_dump<10>("eigenvalues = ", m_eigenvalues);
    m_logger->data_dump("solutions = ", m_solutions);
  }

public:
  //! Set error value for solution *root*
  void set_error(int root, value_type_abs error) override { m_errors.at(root) = error; }
  void set_error(const std::vector<int>& roots, const std::vector<value_type_abs>& errors) override {
    for (size_t i = 0; i < roots.size(); ++i)
      set_error(roots[i], errors[i]);
  }

  const Matrix<value_type>& solutions() const override { return m_solutions; }
  const std::vector<value_type>& eigenvalues() const override { return m_eigenvalues; }
  const std::vector<std::pair<std::size_t, value_type>>& imag_eigval_components() const { return m_imag_eigval_comps; }
  /*!
   * @brief The other member of the complex conjugate eigenvalue pair that root belongs to, if any
   * @note The partner can be >= size(), if it was truncated away
   */
  std::optional<std::size_t> complex_partner(std::size_t root) const {
    const auto it = m_complex_partners.find(root);
    if (it == m_complex_partners.end())
      return std::nullopt;
    return it->second;
  }
  /*!
   * @brief The solutions, followed by the complex partner of the last root if that partner was truncated away
   *
   * Row r holds the coefficients of root r, so the row of a truncated partner is complex_partner(size() - 1).
   */
  const Matrix<value_type>& retained_solutions() const override {
    return m_retained_solutions.empty() ? m_solutions : m_retained_solutions;
  }
  const std::vector<value_type_abs>& errors() const override { return m_errors; }

  //! Number of solutions
  size_t size() const override { return m_solutions.rows(); }

  void set_logger(std::shared_ptr<Logger> logger) override { m_logger = std::move(logger); }

  // FIXME What difference does it make?
  //! Set Hermiticity of the subspace.
  void set_hermiticity(bool hermitian) { m_hermitian = hermitian; }
  bool get_hermiticity() { return m_hermitian; }
  //! Set value of augmented hessian parameter. If 0, than augmented Hessian is not used.
  void set_augmented_hessian(value_type_abs parameter) { m_augmented_hessian = parameter; }
  value_type_abs get_augmented_hessian() { return m_augmented_hessian; }

protected:
  Matrix<value_type> m_solutions;                                      //!< solution matrix with row vectors
  std::vector<value_type> m_eigenvalues;                               //!< eigenvalues
  std::vector<std::pair<std::size_t, value_type>> m_imag_eigval_comps; //!< eigenvalues
  std::map<std::size_t, std::size_t> m_complex_partners;               //!< root -> other member of its complex pair
  Matrix<value_type> m_retained_solutions; //!< solutions and truncated partner of the last root, empty without one
  std::vector<value_type_abs> m_errors;                                //!< errors in subspace solutions
  std::shared_ptr<Logger> m_logger{};

public:
  //! threshold to select null space during SVD in eigenproblem, rescaled from double to the working precision
  value_type_abs m_svd_solver_threshold = precision_scaled<value_type_abs>(1e-14);
protected:
  bool m_hermitian = false;       //!< flags the matrix as Hermitian
  //! value of augmented hessian parameter. If 0, than augmented Hessian is not used
  value_type_abs m_augmented_hessian = 0;
};

} // namespace molpro::linalg::itsolv::subspace

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SUBSPACESOLVERLINEIG_H
