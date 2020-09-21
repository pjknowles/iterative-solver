#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_QSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_QSPACE_H
#include <cassert>
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
  //! previous parameters).
  //! @Returns scaling constants for this and x leading to result, i.e. a and b in \f$ r' = a r + b x \f$
  std::pair<double, double> merge(const QParam<Q>& x, array::ArrayHandler<Q, Q>& handler) {
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
    return {normalisation_constant, a * normalisation_constant};
  }
};

//! Generates vector of reference wrappers to param and action in a range of QParam objects
template <class Q, class ForwardIt>
auto wrap_params(ForwardIt begin, ForwardIt end) {
  using VecRefQ = std::vector<std::reference_wrapper<Q>>;
  auto param_action = std::array<VecRefQ, 2>{};
  for (auto it = begin; it != end; ++it) {
    param_action[0].emplace_back(*it->param);
    param_action[1].emplace_back(*it->action);
  }
  return param_action;
}

//! Generate new difference vectors based on current and last working set
template <class R, class Q, class P>
std::pair<std::list<QParam<Q>>, std::vector<size_t>>
update(R& qparam, R& qaction, const std::vector<std::reference_wrapper<R>>& params,
       const std::vector<std::reference_wrapper<R>>& actions, const std::vector<Q>& last_params,
       const std::vector<Q>& last_actions, const std::vector<size_t>& working_set, ArrayHandlers<R, Q, P>& handlers) {
  assert(params.size() == last_params.size() && params.size() == actions.size() &&
         last_params.size() == last_actions.size() && "Must provide consistent number of input parameters");
  auto used_working_set = std::vector<size_t>{};
  auto qparams = std::list<QParam<Q>>{};
  for (size_t i = 0; i < working_set.size(); ++i) {
    handlers.rq().copy(qparam, params.at(i));
    handlers.rq().copy(qaction, actions.at(i));
    handlers.rq().axpy(-1, last_params.at(i), qparam);
    handlers.rq().axpy(-1, last_actions.at(i), qaction);
    auto qq = handlers.rr().dot(qparam, qparam);
    auto norm = std::sqrt(qq);
    // FIXME no orthogonalisation is done
    // FIXME what do we do if norm is very small?
    if (norm > 1.0e-14) {
      handlers.rr().scal(1. / norm, qparam);
      handlers.rr().scal(1. / norm, qaction);
      auto&& q = qspace::QParam<Q>{std::make_unique<Q>(handlers.qr().copy(qparam)),
                                   std::make_unique<Q>(handlers.qr().copy(qaction)),
                                   working_set[i],
                                   false,
                                   1. / norm,
                                   1};
      qparams.emplace_back(std::move(q));
      used_working_set.emplace_back(working_set[i]);
    }
  }
  return {std::move(qparams), used_working_set};
}

//! Updates equation data in the QxQ part of the subspace
template <class Q>
void update_qq_subspace(const std::vector<std::reference_wrapper<Q>>& old_params,
                        const std::vector<std::reference_wrapper<Q>>& old_actions,
                        const std::vector<std::reference_wrapper<Q>>& new_params,
                        const std::vector<std::reference_wrapper<Q>>& new_actions, SubspaceData& data,
                        array::ArrayHandler<Q, Q>& handler) {
  auto nQold = old_params.size();
  auto nQnew = new_params.size();
  auto nQ = nQnew + nQold;
  data[EqnData::S].resize({nQ, nQ});
  data[EqnData::H].resize({nQ, nQ});
  auto ov_params_old_new = util::overlap(old_params, new_params, handler);
  auto ov_actions_old_new = util::overlap(old_actions, new_actions, handler);
  auto ov_params_new_old = util::overlap(new_params, old_params, handler);
  auto ov_actions_new_old = util::overlap(new_actions, old_actions, handler);
  auto ov_params_new_new = util::overlap(new_params, new_params, handler);
  auto ov_actions_new_new = util::overlap(new_actions, new_actions, handler);
  data[EqnData::S].slice({nQold, 0}, {nQ, nQold}) = ov_params_new_old;
  data[EqnData::H].slice({nQold, 0}, {nQ, nQold}) = ov_actions_new_old;
  data[EqnData::S].slice({0, nQold}, {nQold, nQ}) = ov_params_old_new;
  data[EqnData::H].slice({0, nQold}, {nQold, nQ}) = ov_actions_old_new;
  data[EqnData::S].slice({nQold, nQold}, {nQ, nQ}) = ov_params_new_new;
  data[EqnData::H].slice({nQold, nQold}, {nQ, nQ}) = ov_actions_new_new;
}

//! Updates equation data in the RxQ part of the subspace
template <class R, class Q>
void update_qr_subspace(const std::vector<std::reference_wrapper<Q>>& qparams,
                        const std::vector<std::reference_wrapper<Q>>& qactions,
                        const std::vector<std::reference_wrapper<R>>& rparams,
                        const std::vector<std::reference_wrapper<R>>& ractions, SubspaceData& qr, SubspaceData& rq,
                        array::ArrayHandler<Q, R>& handler_qr, array::ArrayHandler<R, Q>& handler_rq) {
  auto nQ = qparams.size();
  auto nR = rparams.size();
  qr[EqnData::S].resize({nQ, nR});
  qr[EqnData::H].resize({nQ, nR});
  rq[EqnData::S].resize({nR, nQ});
  rq[EqnData::H].resize({nR, nQ});
  auto ov_params_qr = util::overlap(qparams, rparams, handler_qr);
  auto ov_actions_qr = util::overlap(qactions, ractions, handler_qr);
  // FIXME in hermitian cases rq is redundant
  auto ov_params_rq = util::overlap(rparams, qparams, handler_rq);
  auto ov_actions_rq = util::overlap(ractions, qactions, handler_rq);
  qr[EqnData::S].slice() = ov_params_qr;
  qr[EqnData::H].slice() = ov_actions_qr;
  rq[EqnData::S].slice() = ov_params_rq;
  rq[EqnData::H].slice() = ov_actions_rq;
}

//! Merges subspace data for q params i and j \f$ q_i = a q_i + b q_j \f$
void merge_subspace_qq(size_t i, size_t j, double a, double b, SubspaceData& qq);

//! Merges subspace data for q params i and j \f$ q_i = a q_i + b q_j \f$
void merge_subspace_qr(size_t i, size_t j, double a, double b, SubspaceData& qr);

//! Merges subspace data for q params i and j \f$ q_i = a q_i + b q_j \f$
void merge_subspace_rq(size_t i, size_t j, double a, double b, SubspaceData& rq);

//! Removes data associates with q parameter i from qq, qr and rq blocks
void erase_subspace(size_t i, SubspaceData& qq, SubspaceData& qr, SubspaceData& rq);

} // namespace qspace

/*!
 * @brief Container for building and managing the Q space parameters.
 *
 * **Add description of what Q space is**
 *
 * This QSpace is built as difference of current solutions from the previous.
 * The difference vectors are normalised and can be optionally orthogonalised (not implemented yet).
 *
 * **Add documentation on merging or removing subspace elements**
 *
 * @tparam R array for R space
 * @tparam Q array for Q space
 * @tparam P array for P space
 */
template <class R, class Q, class P>
struct QSpace {
  using VecRefR = typename RSpace<R, Q, P>::VecRefR;
  using VecRefQ = std::vector<std::reference_wrapper<Q>>;

  SubspaceData data = null_data<EqnData::H, EqnData::S>(); //!< QxQ block of subspace data
  SubspaceData qr = null_data<EqnData::H, EqnData::S>();   //!< QxR block of subspace data
  SubspaceData rq = null_data<EqnData::H, EqnData::S>();   //!< RxQ block of subspace data

  explicit QSpace(std::shared_ptr<ArrayHandlers<R, Q, P>> handlers) : m_handlers(std::move(handlers)) {}

  void update(const RSpace<R, Q, P>& rs, IterativeSolver<R, Q, P>& solver) {
    auto& dummy = rs.dummy(2);
    auto result = qspace::update(dummy.at(0), dummy.at(1), rs.params(), rs.actions(), rs.last_params(),
                                 rs.last_actions(), rs.working_set(), *m_handlers);
    auto& new_qparams = result.first;
    m_used_working_set = result.second;
    auto old_params_actions = qspace::wrap_params<Q>(m_params.begin(), m_params.end());
    auto new_params_actions = qspace::wrap_params<Q>(new_qparams.begin(), new_qparams.end());
    m_params.splice(m_params.end(), new_qparams);
    auto all_params_actions = qspace::wrap_params<Q>(m_params.begin(), m_params.end());
    qspace::update_qq_subspace(old_params_actions[0], old_params_actions[1], new_params_actions[0],
                               new_params_actions[1], data, m_handlers->qq());
    qspace::update_qr_subspace(all_params_actions[0], all_params_actions[1], rs.params(), rs.actions(), qr, rq,
                               m_handlers->rq(), m_handlers->qr());
  }

  void add_converged(const VecRefR& params, VecRefR& actions, const std::vector<size_t>& roots) {
    for (size_t i = 0; i < params.size(); ++i) {
      auto&& q = qspace::QParam<Q>{std::make_unique<Q>(m_handlers.qr().copy(params[i])),
                                   std::make_unique<Q>(m_handlers.qr().copy(actions[i])),
                                   roots[i],
                                   true,
                                   1,
                                   1};
      m_params.emplace_back(std::move(q));
    }
    auto nQ = m_params.size();
    auto nQnew = params.size();
    auto nQprev = nQ - nQnew;
    auto old_params_actions = qspace::wrap_params<Q>(m_params.begin(), next(m_params.begin(), nQprev));
    auto new_params_actions = qspace::wrap_params<Q>(next(m_params.begin(), nQprev), m_params.end());
    qspace::update_qq_subspace(old_params_actions[0], old_params_actions[1], new_params_actions[0],
                               new_params_actions[1], data, m_handlers->qq());
  }

  //! Vector of root indices for r vectors that were used to generate new q vectors. Converged solutions are not
  //! included.
  auto& used_working_set() const { return m_used_working_set; }

  //! Returns indices of q parameters corresponding to root that can be modified. Converged solutions and latest q
  //! vector for that root are excluded.
  auto modification_candidates(size_t root) const {
    auto candidates = std::vector<size_t>{};
    size_t i = 0;
    for (auto it = m_params.begin(); it != m_params.end(); ++it, ++i)
      if (it->root == root)
        candidates.emplace_back(i);
    if (!candidates.empty())
      candidates.resize(candidates.size() - 1);
    return candidates;
  }

  //! Merges a pair of q vectors if they belong to the same root and updates the subspace data
  void merge(const std::pair<size_t, size_t>& pair) {
    assert(m_params.size() > pair.first && m_params.size() > pair.second);
    auto i = std::min(pair.first, pair.second);
    auto j = std::max(pair.first, pair.second);
    if (i == j)
      throw std::runtime_error("attempting to merge the same vector");
    auto left = std::next(begin(m_params), i);
    auto right = std::next(begin(m_params), j);
    auto root = left->root;
    if (root != right->root)
      throw std::runtime_error("attempting to merge difference vectors corresponding to different roots");
    auto first_difference_vector =
        std::find_if(begin(m_params), end(m_params), [root](const auto& el) { return el.root == root; });
    if (left == first_difference_vector) {
      erase(i);
    } else {
      double a, b;
      std::tie(a, b) = left->merge(*right, m_handlers->qq());
      m_params.erase(right);
      qspace::merge_subspace_qq(i, j, a, b, data);
      qspace::merge_subspace_qr(i, j, a, b, qr);
      qspace::merge_subspace_rq(i, j, a, b, rq);
    }
  }

  //! Erases q parameter i. @param i index in the current Q space
  void erase(size_t i) {
    assert(m_params.size() > i);
    m_params.erase(std::next(begin(m_params), i));
    // remove associated rows and columns from the data
    qspace::erase_subspace(i, data, qr, rq);
  }

protected:
  std::shared_ptr<ArrayHandlers<R, Q, P>> m_handlers;
  std::vector<size_t> m_used_working_set; //!< root indices of r vectors that were used to generate new q vectors
  std::list<qspace::QParam<Q>> m_params;  //!< q vectors constructed as differences
  bool m_orthogonalise = true;            //!< whether to orthogonalise a new Q vector relative to its R vector
};

} // namespace subspace
} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_QSPACE_H
