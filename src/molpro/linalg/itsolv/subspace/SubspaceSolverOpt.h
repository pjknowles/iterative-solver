#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SUBSPACESOLVEROPT_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SUBSPACESOLVEROPT_H
#include <molpro/linalg/itsolv/subspace/ISubspaceSolver.h>
#include <molpro/linalg/itsolv/subspace/IXSpace.h>
#include <molpro/linalg/itsolv/subspace/Matrix.h>

namespace molpro::linalg::itsolv::subspace {
template <class T>
std::string as_string(const std::vector<T>& m, int precision = 6) {
  return as_string(Matrix<T>{m, {1, m.size()}});
}


/*!
 * @brief Solves subspace problem for linear eigenvalues and system of linear equations
 */
template <class RT, class QT, class PT>
class SubspaceSolverOpt : public ISubspaceSolver<RT, QT, PT> {
public:
  using value_type = typename ISubspaceSolver<RT, QT, PT>::value_type;
  using value_type_abs = typename ISubspaceSolver<RT, QT, PT>::value_type_abs;
  using R = typename ISubspaceSolver<RT, QT, PT>::R;
  using Q = typename ISubspaceSolver<RT, QT, PT>::Q;
  using P = typename ISubspaceSolver<RT, QT, PT>::P;

  explicit SubspaceSolverOpt(std::shared_ptr<Logger> logger) : m_logger(std::move(logger)) {}

  void solve(IXSpace<R, Q, P>& xspace, const size_t nroots_max) override {
    m_logger->msg("SubspaceSolverOpt::solve", Logger::Trace);
    assert(xspace.data.end() != xspace.data.find(EqnData::value));
    auto values = xspace.data[EqnData::value];
    std::cout << "xspace.size() "<<xspace.size()<<std::endl;
    std::cout << "values.size() "<<values.size()<<std::endl;
    assert(xspace.size() == values.size());

    if (true) {
      solve_steepest(xspace);
    } else {
    }
  }

protected:
  void solve_steepest(IXSpace<R, Q, P>& xspace) {
    m_logger->msg("SubspaceSolverOpt::solve_steepest", Logger::Trace);
    auto h = xspace.data[EqnData::H];
    auto s = xspace.data[EqnData::S];
    auto value = xspace.data[EqnData::value];
    if (m_logger->data_dump) {
      m_logger->msg("S = " + as_string(s), Logger::Info);
      m_logger->msg("H = " + as_string(h, 15), Logger::Info);
      m_logger->msg("value = " + as_string(value, 15), Logger::Info);
    }
    auto dim = h.rows();
    auto evec = std::vector<value_type>{};
    int verbosity = m_logger->max_trace_level == Logger::Info ? 3 : 0;
    m_solutions.resize({1, dim});
    m_solutions.slice().fill(0);
    m_solutions(0,0)=1;
    m_errors.assign(1,h(0,0)); // FIXME
    if (m_logger->data_dump) {
      m_logger->msg("solution = " + as_string(m_solutions), Logger::Info);
    }
  }

public:
  //! Set error value for solution *root*
  void set_error(int root, value_type_abs error) override { m_errors.at(root) = error; }
  void set_error(const std::vector<int>& roots, const std::vector<value_type_abs>& errors) override {
    for (size_t i = 0; i < roots.size(); ++i)
      set_error(roots[i], errors[i]);
  }

  const Matrix<value_type>& solutions() const override { return m_solutions; }
  const std::vector<value_type>& eigenvalues() const override { throw std::logic_error("eigenvalues() not available in non-linear method"); }
  const std::vector<value_type_abs>& errors() const override { return m_errors; }

  //! Number of solutions
  size_t size() const override { return m_solutions.rows(); }


protected:
  Matrix<value_type> m_solutions;        //!< solution matrix with row vectors
  std::vector<value_type_abs> m_errors;  //!< errors in subspace solutions
  std::shared_ptr<Logger> m_logger{};

public:
  value_type_abs m_svd_solver_threshold = 1.0e-14; //!< threshold to select null space during SVD in eigenproblem
};

} // namespace molpro::linalg::itsolv::subspace

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SUBSPACESOLVEROPT_H
