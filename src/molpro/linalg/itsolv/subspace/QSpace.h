#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_QSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_QSPACE_H
#include <cassert>
#include <list>

#include <molpro/linalg/itsolv/subspace/CSpace.h>
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
  std::unique_ptr<Q> param;  //!< parameter
  std::unique_ptr<Q> action; //!< corresponding action
  size_t id;                 //!< unique id
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

//! Prepends equation data for new parameters to the Q subspace data
template <class Q>
void update_qq_subspace(const std::vector<std::reference_wrapper<Q>>& new_params,
                        const std::vector<std::reference_wrapper<Q>>& new_actions,
                        const std::vector<std::reference_wrapper<Q>>& old_params,
                        const std::vector<std::reference_wrapper<Q>>& old_actions, const SubspaceData& qq_new,
                        SubspaceData& data, array::ArrayHandler<Q, Q>& handler) {
  auto nQnew = new_params.size();
  auto nQold = old_params.size();
  auto oQnew = 0;
  auto oQold = oQnew;
  auto nQ = nQnew + nQold;
  auto s = data[EqnData::S];
  auto h = data[EqnData::H];
  s.resize({nQ, nQ});
  h.resize({nQ, nQ});
  s.slice({oQold, oQold}, {oQold + nQold, oQold + nQold}) = data[EqnData::S];
  s.slice({oQold, oQold}, {oQold + nQold, oQold + nQold}) = data[EqnData::H];
  s.slice({oQnew, oQnew}, {oQnew + nQnew, oQnew + nQnew}) = qq_new.at(EqnData::S);
  s.slice({oQnew, oQnew}, {oQnew + nQnew, oQnew + nQnew}) = qq_new.at(EqnData::H);
  s.slice({oQold, oQnew}, {oQold + nQold, oQnew + nQnew}) = util::overlap(old_params, new_params, handler);
  h.slice({oQold, oQnew}, {oQold + nQold, oQnew + nQnew}) = util::overlap(old_params, new_actions, handler);
  s.slice({oQnew, oQold}, {oQnew + nQnew, oQold + nQold}) = util::overlap(new_params, old_params, handler);
  h.slice({oQnew, oQold}, {oQnew + nQnew, oQold + nQold}) = util::overlap(new_params, old_actions, handler);
  data[EqnData::S] = std::move(s);
  data[EqnData::H] = std::move(h);
}

//! Updates equation data in the RxQ part of the subspace
template <class R, class Q>
void update_subspace(const std::vector<std::reference_wrapper<Q>>& qparams,
                     const std::vector<std::reference_wrapper<Q>>& qactions,
                     const std::vector<std::reference_wrapper<R>>& rparams,
                     const std::vector<std::reference_wrapper<R>>& ractions, SubspaceData& qr, SubspaceData& rq,
                     array::ArrayHandler<Q, std::decay_t<R>>& handler_qr,
                     array::ArrayHandler<std::decay_t<R>, Q>& handler_rq) {
  auto nQ = qparams.size();
  auto nR = rparams.size();
  qr[EqnData::S].resize({nQ, nR});
  qr[EqnData::H].resize({nQ, nR});
  rq[EqnData::S].resize({nR, nQ});
  rq[EqnData::H].resize({nR, nQ});
  qr[EqnData::S].slice() = util::overlap(qparams, rparams, handler_qr);
  qr[EqnData::H].slice() = util::overlap(qparams, ractions, handler_qr);
  rq[EqnData::S].slice() = util::overlap(rparams, qparams, handler_rq);
  rq[EqnData::H].slice() = util::overlap(rparams, qactions, handler_rq);
}

//! Removes data associates with q parameter i from qq, qr and rq blocks
void erase_subspace(size_t i, SubspaceData& qq, SubspaceData& qr, SubspaceData& rq, SubspaceData& qc,
                    SubspaceData& cq) {
  for (auto d : {EqnData::S, EqnData::H}) {
    qq[d].remove_row_col(i, i);
    qr[d].remove_row(i);
    rq[d].remove_col(i);
    qc[d].remove_row(i);
    cq[d].remove_col(i);
  }
}

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
template <class Rt, class Qt, class Pt>
struct QSpace {
  using R = Rt;
  using Q = Qt;
  using P = Pt;
  using VecRefQ = std::vector<std::reference_wrapper<Q>>;
  using CVecRefQ = std::vector<std::reference_wrapper<Q>>;

  SubspaceData data = null_data<EqnData::H, EqnData::S>(); //!< QxQ block of subspace data
  SubspaceData qr = null_data<EqnData::H, EqnData::S>();   //!< QxR block of subspace data
  SubspaceData rq = null_data<EqnData::H, EqnData::S>();   //!< RxQ block of subspace data
  SubspaceData qc = null_data<EqnData::H, EqnData::S>();   //!< QxC block of subspace data
  SubspaceData cq = null_data<EqnData::H, EqnData::S>();   //!< CxQ block of subspace data

  explicit QSpace(std::shared_ptr<ArrayHandlers<R, Q, P>> handlers, std::shared_ptr<Logger> logger)
      : m_handlers(std::move(handlers)), m_logger(std::move(logger)) {}

  // FIXME Need to decide whether to remove rs entirely
  // FIXME Need to calculate subspace data using R vectors
  void update(RSpace<R, Q, P>& rs, const CSpace<R, Q, P, double>& ss, IterativeSolver<R, Q, P>& solver) {
    m_logger->msg("QSpace::update", Logger::Trace);
    auto new_qparams = std::list<qspace::QParam<Q>>{};
    for (size_t i = 0; i < rs.size(); ++i) {
      new_qparams.emplace_back(qspace::QParam<Q>{std::make_unique<Q>(m_handlers->qr().copy(rs.params().at(i))),
                                                 std::make_unique<Q>(m_handlers->qr().copy(rs.actions().at(i))),
                                                 m_unique_id++});
    }
    auto old_params_actions = qspace::wrap_params<Q>(m_params.begin(), m_params.end());
    auto new_params_actions = qspace::wrap_params<Q>(new_qparams.begin(), new_qparams.end());
    m_params.splice(m_params.begin(), new_qparams);
    auto all_params_actions = qspace::wrap_params<Q>(m_params.begin(), m_params.end());
    qspace::update_qq_subspace(new_params_actions[0], new_params_actions[1], old_params_actions[0],
                               old_params_actions[1], rs.data, data, m_handlers->qq());
    rs.clear();
    qspace::update_subspace(all_params_actions[0], all_params_actions[1], rs.params(), rs.actions(), qr, rq,
                            m_handlers->rq(), m_handlers->qr());
    qspace::update_subspace(all_params_actions[0], all_params_actions[1], ss.params(), ss.actions(), qc, cq,
                            m_handlers->qq(), m_handlers->qq());
    if (m_logger->data_dump) {
      m_logger->msg("Sqq = " + as_string(data[EqnData::S]), Logger::Info);
      m_logger->msg("Hqq = " + as_string(data[EqnData::H]), Logger::Info);
      m_logger->msg("Sqr = " + as_string(qr[EqnData::S]), Logger::Info);
      m_logger->msg("Hqr = " + as_string(qr[EqnData::H]), Logger::Info);
      m_logger->msg("Srq = " + as_string(rq[EqnData::S]), Logger::Info);
      m_logger->msg("Hrq = " + as_string(rq[EqnData::H]), Logger::Info);
      m_logger->msg("Sqc = " + as_string(qc[EqnData::S]), Logger::Info);
      m_logger->msg("Hqc = " + as_string(qc[EqnData::H]), Logger::Info);
      m_logger->msg("Scq = " + as_string(cq[EqnData::S]), Logger::Info);
      m_logger->msg("Hcq = " + as_string(cq[EqnData::H]), Logger::Info);
    }
  }

  void clear() {
    m_params.clear();
    for (auto d : {EqnData::H, EqnData::S}) {
      data[d].clear();
      qr[d].clear();
      rq[d].clear();
      qc[d].clear();
      cq[d].clear();
    }
  }

  //! Erases q parameter i. @param i index in the current Q space
  void erase(size_t i) {
    assert(m_params.size() > i);
    auto it = std::next(begin(m_params), i);
    m_params.erase(it);
    qspace::erase_subspace(i, data, qr, rq, qc, cq);
  }

  size_t size() const { return m_params.size(); }

  VecRefQ params() {
    auto qparams = VecRefQ{};
    std::transform(begin(m_params), end(m_params), std::back_inserter(qparams),
                   [](auto& q) { return std::ref(*q.param); });
    return qparams;
  }

  CVecRefQ params() const {
    auto qparams = CVecRefQ{};
    std::transform(begin(m_params), end(m_params), std::back_inserter(qparams),
                   [](auto& q) { return std::ref(*q.param); });
    return qparams;
  }

  VecRefQ actions() {
    auto qactions = VecRefQ{};
    std::transform(begin(m_params), end(m_params), std::back_inserter(qactions),
                   [](auto& q) { return std::ref(*q.action); });
    return qactions;
  }

  CVecRefQ actions() const {
    auto qactions = CVecRefQ{};
    std::transform(begin(m_params), end(m_params), std::back_inserter(qactions),
                   [](auto& q) { return std::ref(*q.action); });
    return qactions;
  }

protected:
  std::shared_ptr<ArrayHandlers<R, Q, P>> m_handlers;
  std::shared_ptr<Logger> m_logger;
  size_t m_unique_id{0};                 //!< unique id for any new parameter set
  std::list<qspace::QParam<Q>> m_params; //!< q parameter sets with new parameters first
};

} // namespace subspace
} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_QSPACE_H
