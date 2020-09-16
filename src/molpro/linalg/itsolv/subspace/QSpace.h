#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_QSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_QSPACE_H
#include <list>
#include <molpro/linalg/itsolv/subspace/PSpace.h>
#include <molpro/linalg/itsolv/subspace/RSpace.h>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace subspace {

namespace qspace {
/*!
 * @brief Parameter in the Q space. Can be a difference vector or a converged solution.
 *
 * Q space is constructed by appending differences with previous solution, normalised and orthogonalised (w.r.t working
 * param), and by appending converged solutions.
 *
 * \f$ q = \alpha (r - \beta d)\f$
 *
 * @note I will assume difference vectors are not orthogonalised to make the maths simpler. I'll work it out later.
 */
template <class Q>
struct QParam {
  std::unique_ptr<Q> param;          //!< parameter
  std::unique_ptr<Q> action;         //!< corresponding action
  size_t root = 0;                   //!< corresponding root index
  bool converged = false;            //!< whether this is a converged vector
  double normalisation_constant = 1; //!< scaling constant which makes difference vector normalised, see \f$\alpha\f$
  double orthogonalisation_constant =
      1; //!< scaling constant which makes difference vector orthogonal to r, see \f$\beta\f$

  //! Merge this difference parameter with x, in a way that preserves the difference property (ability to reconstruct
  //! previous parameters)
  void merge(const QParam<Q>& x, array::ArrayHandler<Q, Q>& handler) {
    if (converged || x.converged)
      throw std::runtime_error("attempting to merge converged solutions");
    if (root != x.root)
      throw std::runtime_error("attempting to merge q vectors from different roots");
    auto a = normalisation_constant / x.normalisation_constant;
    handler.axpy(a, *x.param, *param);
    handler.axpy(a, *x.action, *action);
    auto dot = handler.dot(*param, *action);
    normalisation_constant = 1. / std::sqrt(dot);
    handler.scal(normalisation_constant, *param);
    handler.scal(normalisation_constant, *action);
  }
};
} // namespace qspace

template <class R, class Q, class P>
struct QSpace {
  using typename RSpace<R, Q, P>::VecRefR;
  SubspaceData data = null_data<EqnData::H, EqnData::S>();

  void update(const RSpace<R, Q, P>& rs, IterativeSolver<R, Q, P>& solver) {
    m_qr = null_data<EqnData::H, EqnData::S>();
    m_rq = null_data<EqnData::H, EqnData::S>();
    auto& dummy = rs.dummy(2);
    auto& qparam = dummy[0];
    auto& qaction = dummy[1];
    for (size_t i = 0; i < rs.working_set().size(); ++i) {
      auto& param = rs.params()[i];
      auto& action = rs.actions()[i];
      auto& last_param = rs.last_params()[i];
      auto& last_action = rs.last_actions()[i];
      m_handlers->rq().copy(qparam, param);
      m_handlers->rq().axpy(-1, last_param, qparam);
      m_handlers->rq().axpy(-1, last_action, qaction);
      auto qq = m_handlers->rr().dot(qparam, qparam);
      auto norm = 1. / std::sqrt(qq);
      // FIXME no orthogonalisation is done
      // FIXME what do we do if norm is very small?
      if (norm > 1.0e-14) {
        m_handlers->rr().scal(norm, qparam);
        m_handlers->rr().scal(norm, qaction);
        auto&& q = qspace::QParam<Q>{std::make_unique<Q>(m_handlers.qr().copy(qparam)),
                                     std::make_unique<Q>(m_handlers.qr().copy(qaction)),
                                     rs.working_set()[i],
                                     false,
                                     norm,
                                     1};
        m_params.emplace_back(std::move(q));
      }
    }
    // update subspace data. This has to be done with all new parameters for efficiency.
  }

  void add_converged(const R& param, const R& action, size_t root) {
    auto&& q = qspace::QParam<Q>{std::make_unique<Q>(m_handlers.qr().copy(param)),
                                 std::make_unique<Q>(m_handlers.qr().copy(action)),
                                 root,
                                 true,
                                 1,
                                 1};
    m_params.push_back(std::move(q));
    // update qq
  }

protected:
  std::shared_ptr<ArrayHandlers<R, Q, P>> m_handlers;
  SubspaceData m_qr = null_data<EqnData::H, EqnData::S>(); //!< QxR section of subspace data
  SubspaceData m_rq = null_data<EqnData::H, EqnData::S>(); //!< RxQ section of subspace data
  std::list<qspace::QParam<Q>> m_params;                   //!< q vectors constructed as differences
  bool m_orthogonalise = true; //!< whether to orthogonalise a new Q vector relative to its R vector
};

} // namespace subspace
} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_QSPACE_H
