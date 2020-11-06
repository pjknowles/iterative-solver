#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_QSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_QSPACE_H
#include <cassert>
#include <list>
#include <array>

#include <molpro/linalg/itsolv/ArrayHandlers.h>
#include <molpro/linalg/itsolv/Logger.h>
#include <molpro/linalg/itsolv/subspace/Dimensions.h>
#include <molpro/linalg/itsolv/subspace/PSpace.h>
#include <molpro/linalg/itsolv/wrap.h>

namespace molpro::linalg::itsolv::subspace {

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

  /*!
   * @brief Prepends parameters to the start of Q space
   * @param params new parameters
   * @param actions new actions
   * @param qq Equation data block between new Q parameters
   * @param qx data block between new Q parameters and current X space
   * @param xq data block between current X space and the new Q parameters
   * @param dims current dimensions
   * @param old_data current data
   */
  void update(const CVecRef<R>& params, const CVecRef<R>& actions, const SubspaceData& qq, const SubspaceData& qx,
              const SubspaceData& xq, const xspace::Dimensions& dims, SubspaceData& old_data) {
    m_logger->msg("QSpace::update", Logger::Trace);
    auto it_begin = m_params.begin();
    for (size_t i = 0; i < params.size(); ++i) {
      m_params.emplace(it_begin,
                       qspace::QParam<Q>{std::make_unique<Q>(m_handlers->qr().copy(params.at(i))),
                                         std::make_unique<Q>(m_handlers->qr().copy(actions.at(i))), m_unique_id++});
    }
    const size_t nQnew = params.size();
    const auto nXnew = dims.nX + nQnew;
    auto data = old_data;
    for (auto d : {EqnData::H, EqnData::S}) {
      data[d].resize({dims.nX + nQnew, dims.nX + nQnew});
      data[d].slice({dims.oQ + nQnew, dims.oQ + nQnew}, {data[d].rows(), data[d].cols()}) =
          old_data[d].slice({dims.oQ, dims.oQ}, {dims.nX, dims.nX});
      data[d].slice({dims.oQ, dims.oQ}, {dims.oQ + nQnew, dims.oQ + nQnew}) = qq.at(d);
      data[d].slice({dims.oQ, 0}, {dims.oQ + nQnew, dims.oQ}) = qx.at(d).slice({0, 0}, {nQnew, dims.oQ});
      data[d].slice({dims.oQ, dims.oQ + nQnew}, {dims.oQ + nQnew, nXnew}) =
          qx.at(d).slice({0, dims.oQ}, {nQnew, dims.nX});
      data[d].slice({0, dims.oQ}, {dims.oQ, dims.oQ + nQnew}) = xq.at(d).slice({0, 0}, {dims.oQ, nQnew});
      data[d].slice({dims.oQ + nQnew, dims.oQ}, {nXnew, dims.oQ + nQnew}) =
          xq.at(d).slice({dims.oQ, 0}, {dims.nX, nQnew});
      old_data[d] = data[d];
    }
    if (m_logger->data_dump) {
      m_logger->msg("S = " + as_string(data.at(EqnData::S)), Logger::Info);
      m_logger->msg("H = " + as_string(data.at(EqnData::H)), Logger::Info);
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

} // namespace molpro::linalg::itsolv::subspace

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_QSPACE_H
