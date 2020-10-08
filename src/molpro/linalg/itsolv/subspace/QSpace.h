#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_QSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_QSPACE_H
#include <cassert>
#include <list>

#include <molpro/linalg/itsolv/subspace/CSpace.h>
#include <molpro/linalg/itsolv/subspace/PSpace.h>
#include <molpro/linalg/itsolv/subspace/RSpace.h>
#include <molpro/linalg/itsolv/wrap.h>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace subspace {

namespace qspace {
/*!
 * @brief Parameter in the Q space.
 */
template <class Q>
struct QParam {
  std::unique_ptr<Q> param;  //!< parameter
  std::unique_ptr<Q> action; //!< corresponding action
  size_t id;                 //!< unique id
};

//! Generates vector of reference wrappers to param and action in a range of QParam objects
template <class Q, class ForwardIt>
auto cwrap_params(ForwardIt begin, ForwardIt end) {
  auto param_action = std::array<CVecRef<Q>, 2>{};
  for (auto it = begin; it != end; ++it) {
    param_action[0].emplace_back(*it->param);
    param_action[1].emplace_back(*it->action);
  }
  return param_action;
}

template <class Q, class ForwardIt>
auto wrap_params(ForwardIt begin, ForwardIt end) {
  auto param_action = std::array<VecRef<Q>, 2>{};
  for (auto it = begin; it != end; ++it) {
    param_action[0].emplace_back(*it->param);
    param_action[1].emplace_back(*it->action);
  }
  return param_action;
}

//! Prepends equation data for new parameters to the Q subspace data
template <class Q>
void update_qq_subspace(const CVecRef<Q>& new_params, const CVecRef<Q>& new_actions, const CVecRef<Q>& old_params,
                        const CVecRef<Q>& old_actions, const SubspaceData& qq_new, SubspaceData& data,
                        array::ArrayHandler<Q, Q>& handler) {
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
  h.slice({oQold, oQold}, {oQold + nQold, oQold + nQold}) = data[EqnData::H];
  s.slice({oQnew, oQnew}, {oQnew + nQnew, oQnew + nQnew}) = qq_new.at(EqnData::S);
  h.slice({oQnew, oQnew}, {oQnew + nQnew, oQnew + nQnew}) = qq_new.at(EqnData::H);
  s.slice({oQold, oQnew}, {oQold + nQold, oQnew + nQnew}) = util::overlap(old_params, new_params, handler);
  h.slice({oQold, oQnew}, {oQold + nQold, oQnew + nQnew}) = util::overlap(old_params, new_actions, handler);
  s.slice({oQnew, oQold}, {oQnew + nQnew, oQold + nQold}) = util::overlap(new_params, old_params, handler);
  h.slice({oQnew, oQold}, {oQnew + nQnew, oQold + nQold}) = util::overlap(new_params, old_actions, handler);
  data[EqnData::S] = std::move(s);
  data[EqnData::H] = std::move(h);
}

//! Updates equation data in the RxQ part of the subspace
template <class R, class Q>
void update_subspace(const CVecRef<Q>& qparams, const CVecRef<Q>& qactions, const CVecRef<R>& rparams,
                     const CVecRef<R>& ractions, SubspaceData& qr, SubspaceData& rq,
                     array::ArrayHandler<Q, R>& handler_qr, array::ArrayHandler<R, Q>& handler_rq) {
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

//! Updates qc and cq blocks in CSpace
template <class Q>
void update_cspace(const CVecRef<Q>& params, const CVecRef<Q>& actions, SubspaceData& qc, SubspaceData& cq,
                   const CVecRef<Q>& cparams, const CVecRef<Q>& cactions, array::ArrayHandler<Q, Q>& handler) {
  const auto nQ_new = params.size();
  const auto nC = cq[EqnData::S].rows();
  auto cq_new = cq;
  auto qc_new = qc;
  for (auto d : {EqnData::S, EqnData::H}) {
    cq_new[d].resize({cq[d].rows(), cq[d].cols() + nQ_new});
    cq_new[d].slice({0, nQ_new}, {cq[d].rows(), cq[d].cols() + nQ_new}) = cq[d];
    qc_new[d].resize({qc[d].rows() + nQ_new, qc[d].cols()});
    qc_new[d].slice({nQ_new, 0}, {qc[d].rows() + nQ_new, qc[d].cols()}) = qc[d];
  }
  auto s_qc = util::overlap(params, cparams, handler);
  auto h_qc = util::overlap(params, cactions, handler);
  auto h_cq = util::overlap(cparams, actions, handler);
  qc_new[EqnData::S].slice({0, 0}, {nQ_new, nC}) = s_qc;
  qc_new[EqnData::H].slice({0, 0}, {nQ_new, nC}) = h_qc;
  cq_new[EqnData::H].slice({0, 0}, {nC, nQ_new}) = h_cq;
  for (size_t i = 0; i < nC; ++i)
    for (size_t j = 0; j < nQ_new; ++j)
      cq_new[EqnData::S](i, j) = qc_new[EqnData::S](j, i);
  for (auto d : {EqnData::S, EqnData::H}) {
    cq[d] = std::move(cq_new[d]);
    qc[d] = std::move(qc_new[d]);
  }
}

//! Removes data associates with q parameter i from qq, qr and rq blocks
void erase_subspace(size_t i, SubspaceData& qq, SubspaceData& qr, SubspaceData& rq) {
  for (auto d : {EqnData::S, EqnData::H}) {
    qq[d].remove_row_col(i, i);
    qr[d].remove_row(i);
    rq[d].remove_col(i);
  }
}

} // namespace qspace

/*!
 * @brief Container storing the Q space parameters.
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
  using value_type = typename array::ArrayHandler<R, R>::value_type;
  using value_type_abs = typename array::ArrayHandler<R, R>::value_type_abs;

  explicit QSpace(std::shared_ptr<ArrayHandlers<R, Q, P>> handlers, std::shared_ptr<Logger> logger)
      : m_handlers(std::move(handlers)), m_logger(std::move(logger)) {}

  //! Prepends parameters to the start of Q space
  void update(const CVecRef<R>& params, const CVecRef<R>& actions) {
    m_logger->msg("QSpace::update", Logger::Trace);
    auto it_begin = m_params.begin();
    for (size_t i = 0; i < params.size(); ++i) {
      m_params.emplace(it_begin,
                       qspace::QParam<Q>{std::make_unique<Q>(m_handlers->qr().copy(params.at(i))),
                                         std::make_unique<Q>(m_handlers->qr().copy(actions.at(i))), m_unique_id++});
    }
  }

  void clear() { m_params.clear(); }

  //! Erases q parameter i. @param i index in the current Q space
  void erase(size_t i) {
    assert(m_params.size() > i);
    auto it = std::next(begin(m_params), i);
    m_params.erase(it);
  }

  size_t size() const { return m_params.size(); }

  VecRef<Q> params() {
    auto qparams = qspace::wrap_params<Q>(m_params.begin(), m_params.end());
    return qparams[0];
  }

  CVecRef<Q> params() const {
    auto qparams = qspace::cwrap_params<Q>(m_params.begin(), m_params.end());
    return qparams[0];
  }

  CVecRef<Q> cparams() const { return params(); }

  VecRef<Q> actions() {
    auto qparams = qspace::wrap_params<Q>(m_params.begin(), m_params.end());
    return qparams[1];
  }

  CVecRef<Q> actions() const {
    auto qparams = qspace::cwrap_params<Q>(m_params.begin(), m_params.end());
    return qparams[1];
  }

  CVecRef<Q> cactions() const { return actions(); }

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
