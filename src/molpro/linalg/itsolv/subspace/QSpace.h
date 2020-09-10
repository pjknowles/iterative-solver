#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_QSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_QSPACE_H
#include <list>
#include <molpro/linalg/itsolv/detail/PSpace.h>
#include <molpro/linalg/itsolv/detail/RSpace.h>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace detail {

namespace qspace {

template <class R, class Q>
std::pair<bool, std::pair<double, double>>
generate_new_q(const std::map<size_t, std::reference_wrapper<R>>& params,
               const std::map<size_t, std::reference_wrapper<R>>& actions, const std::vector<Q>& old_params,
               const std::vector<Q>& old_actions, R& qparam, R& qaction, bool do_orthogonalise) {
  // get diff and scaling constants
  // abort if needed
  // form qvec
}

template <class R, class Q>
void expand_subspace(SubspaceData& qq, SubspaceData& qr, SubspaceData& rq, const R& param, const R& action,
                     const std::vector<Q>& params, const std::vector<Q>& actions) {
  // build qq subspace
  // build qr subspace
  // build rq subspace
}

} // namespace qspace

template <class R, class Q, class P>
struct QSpace {
  using typename RSpace<R, Q>::VecRefR;
  SubspaceData subspace = null_data<EqnData::H, EqnData::S>();

  void update(const RSpace<R, Q>& rs, IterativeSolver<R, Q, P>& solver) {
    for (auto root : solver.working_set()) {
      auto qvecs = rs.dummy();
      update(rs.params().at(root), rs.actions().at(root), rs.last_params().at(root), rs.last_actions().at(root),
             qvecs[0], qvecs[1], m_orthogonalise);
    }
  }

  void add_vector(const R& qparam, const R& qaction, size_t root) {
    qspace::expand_subspace(subspace, m_qr, m_rq, m_qparams, m_qactions, qparam, qaction);
    m_qparams.emplace(m_handlers->qr().copy(qparam));
    m_qactions.emplace(m_handlers->qr().copy(qaction));
    m_roots.emplace_back(root);
  }

protected:
  void update(const R& param, const R& action, const Q& last_param, const Q& last_action, R& qparam, R& qaction,
              size_t root) {
    auto optional_factors =
        qspace::generate_new_q(param, action, last_param, last_action, qparam, qaction, m_orthogonalise);
    if (optional_factors.first) {
      // store scaling factors
      add_vector(qparam, qaction, root);
    }
  }

  std::shared_ptr<ArrayHandlers<R, Q, P>> m_handlers;
  SubspaceData m_qr = null_data<EqnData::H, EqnData::S>(); //!< QxR section of subspace data
  SubspaceData m_rq = null_data<EqnData::H, EqnData::S>(); //!< RxQ section of subspace data
  std::list<Q> m_qparams;                                  //!< Q vectors forming the subspace
  std::list<Q> m_qactions;                                 //!< action corresponding to each Q vector
  std::vector<size_t> m_roots;                             //!< for each q vector stores the corresponding root index
  bool m_orthogonalise = true; //!< whether to orthogonalise a new Q vector relative to its R vector
  /*
   * Q vectors will be added, than merged or deleted.
   * The subspace will need to modify the relevant rows and columns
   * subspace.row(i)
   */
};

template <class R, class Q, class P>
struct QSpaceLE : public QSpace<R, Q, P> {
  using QSpace<R, Q, P>::subspace;
  QSpaceLE() : QSpace<R, Q, P>() { subspace = null_data<EqnData::H, EqnData::S, EqnData::rhs>; }

  void update(const RSpace<R, Q>& rs, LinearEigensystem<R, Q, P>& solver) { QSpace<R, Q, P>::update(rs, solver); }
};

} // namespace detail
} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_QSPACE_H
