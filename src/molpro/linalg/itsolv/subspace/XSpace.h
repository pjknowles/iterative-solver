#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_XSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_XSPACE_H
#include <cassert>
#include <molpro/linalg/itsolv/helper.h>
#include <molpro/linalg/itsolv/subspace/CSpace.h>
#include <molpro/linalg/itsolv/subspace/PSpace.h>
#include <molpro/linalg/itsolv/subspace/QSpace.h>
#include <molpro/linalg/itsolv/subspace/XSpaceI.h>
#include <molpro/linalg/itsolv/subspace/build_subspace.h>
#include <molpro/linalg/itsolv/subspace/check_conditioning.h>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace subspace {

template <class R, class Q, class P>
class XSpace : public XSpaceI<R, Q, P> {
public:
  using typename XSpaceI<R, Q, P>::value_type;
  using typename XSpaceI<R, Q, P>::value_type_abs;

  explicit XSpace(const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers, const std::shared_ptr<Logger>& logger)
      : pspace(), qspace(handlers, logger), cspace(handlers, logger), m_logger(logger) {
    this->data = null_data<EqnData::H, EqnData::S>();
  };

  void update_qspace(const std::vector<R>& params, const std::vector<R>& actions) override {
    qspace.update(wrap(params), wrap(actions));
    // update data
  }

  void update_cspace(const std::vector<unsigned int>& roots, const std::vector<R>& params,
                     const std::vector<R>& actions, const std::vector<value_type>& errors) override {
    cspace.update(roots, params, actions, m_eval, this->data, m_evec, m_dim, errors);
    // update data
  }

  void check_conditioning() override {
    auto nX_on_entry = m_dim.nX;
    m_logger->msg("XSpaceLinEig::check_conditioning size of x space = " + std::to_string(m_dim.nX), Logger::Trace);
    m_logger->msg("size of x space before conditioning = " + std::to_string(m_dim.nX), Logger::Debug);
    if (m_logger->data_dump) {
      m_logger->msg("on entry", Logger::Info);
      m_logger->msg("Sxx = " + as_string(this->data[EqnData::S]), Logger::Info);
      m_logger->msg("Hxx = " + as_string(this->data[EqnData::H]), Logger::Info);
    }
    xspace::check_conditioning_gram_schmidt(*this, m_subspace_transformation, m_norm_stability_threshold, *m_logger);
    m_logger->msg("size of x space after conditioning = " + std::to_string(m_dim.nX), Logger::Debug);
    if (m_logger->data_dump && m_dim.nX != nX_on_entry) {
      m_logger->msg("Sxx = " + as_string(this->data[EqnData::S]), Logger::Info);
      m_logger->msg("Hxx = " + as_string(this->data[EqnData::H]), Logger::Info);
    }
  }

  //! Solves the underlying problem in the subspace. solver can be used to pass extra parameters defining the problem
  void solve(const IterativeSolver<R, Q, P>& solver) {
    assert("XSpaceLinEig can only be used with LinearEigensystem solver");
  };

  //! Solves the linear eigensystem problem in the subspace.
  void solve(const LinearEigensystem<R, Q, P>& solver) {
    m_logger->msg("XSpaceLinEig::solve", Logger::Trace);
    check_conditioning();
    auto& h = this->data[EqnData::H];
    auto& s = this->data[EqnData::S];
    if (m_hermitian)
      util::matrix_symmetrize(h);
    auto dim = h.rows();
    auto evec = std::vector<value_type>{};
    itsolv::eigenproblem(evec, m_eval, h.data(), s.data(), dim, m_hermitian, m_svd_solver_threshold, 0);
    auto n_solutions = evec.size() / dim;
    auto full_matrix = Matrix<value_type>{std::move(evec), {n_solutions, dim}};
    auto nroots = std::min(solver.n_roots(), n_solutions);
    m_eval.resize(nroots);
    m_evec.resize({nroots, dim});
    m_evec.slice() = full_matrix.slice({0, 0}, {nroots, dim});
    if (m_logger->data_dump) {
      m_logger->msg("eigenvalues = ", begin(m_eval), end(m_eval), Logger::Debug);
      m_logger->msg("eigenvectors = " + as_string(m_evec), Logger::Info);
    }
  }

  const std::vector<value_type>& eigenvalues() const override { return m_eval; };

  const Matrix<value_type>& solutions() const override { return m_evec; };

  const xspace::Dimensions& dimensions() const override { return m_dim; }

  void erase(size_t i) override {
    if (m_dim.oP >= i && i < m_dim.oP + m_dim.nP) {
      erasep(i - m_dim.oP);
    } else if (m_dim.oQ >= i && i < m_dim.oQ + m_dim.nQ) {
      eraseq(i - m_dim.oQ);
    } else if (m_dim.oC >= i && i < m_dim.oC + m_dim.nC) {
      erasec(i - m_dim.oC);
    }
  }

  void eraseq(size_t i) override {
    qspace.erase(i);
    // remove data
  }

  void erasep(size_t i) override {
    pspace.erase(i);
    // remove data
  }

  void erasec(size_t i) override {
    cspace.erase(i);
    // remove data
  }

  VecRef<P> paramsp() override { return pspace.params(); }
  VecRef<P> actionsp() override { return pspace.actions(); }
  VecRef<Q> paramsq() override { return qspace.params(); }
  VecRef<Q> actionsq() override { return qspace.actions(); }
  VecRef<Q> paramsc() override { return cspace.params(); }
  VecRef<Q> actionsc() override { return cspace.actions(); }

  CVecRef<P> paramsp() const override { return pspace.params(); }
  CVecRef<P> actionsp() const override { return pspace.actions(); }
  CVecRef<Q> paramsq() const override { return qspace.params(); }
  CVecRef<Q> actionsq() const override { return qspace.actions(); }
  CVecRef<Q> paramsc() const override { return cspace.cparams(); }
  CVecRef<Q> actionsc() const override { return cspace.cactions(); }

  CVecRef<P> cparamsp() const override { return pspace.cparams(); }
  CVecRef<P> cactionsp() const override { return pspace.cactions(); }
  CVecRef<Q> cparamsq() const override { return qspace.cparams(); }
  CVecRef<Q> cactionsq() const override { return qspace.cactions(); }
  CVecRef<Q> cparamsc() const override { return cspace.cparams(); }
  CVecRef<Q> cactionsc() const override { return cspace.cactions(); }

  PSpace<R, P> pspace;
  QSpace<R, Q, P> qspace;
  CSpace<R, Q, P, value_type> cspace;

protected:
  std::shared_ptr<Logger> m_logger;
  xspace::Dimensions m_dim;
  bool m_hermitian = false; //!< whether the matrix is Hermitian
  double m_norm_stability_threshold =
      1.0e-5; //!< norm subspace vector after orhtogonalisation must be less than this to trigger removal
  double m_svd_solver_threshold = 1.0e-14; //!< threshold to remove the null space during solution
  Matrix<value_type>
      m_subspace_transformation;  //!< linear transformation of subspace vectors that leads to stable overlap
  Matrix<value_type> m_evec;      //!< eigenvectors stored as columns with ascending eigenvalue
  std::vector<value_type> m_eval; //!< eigenvalues in ascending order
};

} // namespace subspace
} // namespace itsolv
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_XSPACE_H
