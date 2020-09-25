#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_RSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_RSPACE_H
#include <functional>
#include <memory>
#include <molpro/linalg/itsolv/ArrayHandlers.h>
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/Logger.h>
#include <molpro/linalg/itsolv/subspace/SubspaceData.h>
#include <molpro/linalg/itsolv/subspace/util.h>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace subspace {

namespace rspace {
//! Assigns new parameters to previous based on maximum overlap. Ordering is by increasing root indices.
template <class R, class Q>
std::vector<size_t> assign_last_parameters_to_new(const std::vector<Q>& last_params, const std::vector<R>& new_params,
                                                  array::ArrayHandler<Q, R>& handler) {
  using util::overlap;
  using util::wrap;
  assert(new_params.size() >= last_params.size());
  auto ov = overlap(wrap(last_params), wrap(new_params), handler);
  auto row_to_last = std::vector<size_t>(last_params.size());
  std::iota(begin(row_to_last), end(row_to_last), size_t{0});
  auto last_to_new = std::vector<size_t>(last_params.size());
  while (!ov.empty()) {
    auto imax = std::max_element(begin(ov.data()), end(ov.data()));
    auto row_col = ov.to_coord(distance(ov.data().begin(), imax));
    auto last_ind = row_to_last[row_col.first];
    auto new_ind = row_col.second;
    last_to_new[last_ind] = new_ind;
    ov.remove_row(row_col.first);
    row_to_last.erase(row_to_last.begin() + row_col.first);
  }
  return last_to_new;
}

} // namespace rspace

//! Space for the working set of vectors
template <class Rt, class Qt, class Pt>
class RSpace {
public:
  using R = Rt;
  using Q = Qt;
  using P = Pt;
  using VecRefR = std::vector<std::reference_wrapper<R>>;

  //! Matrix and overlap data mapped to the subspace
  SubspaceData data = null_data<EqnData::H, EqnData::S>();

  explicit RSpace(std::shared_ptr<ArrayHandlers<R, Q, P>> handlers, std::shared_ptr<Logger> logger)
      : m_handlers(std::move(handlers)), m_logger(std::move(logger)) {}

  void update(std::vector<R>& parameters, std::vector<R>& actions, IterativeSolver<R, Q, P>& solver) {
    auto ind_last_param_to_new = std::vector<size_t>{};
    if (m_last_params.empty()) {
      m_logger->msg("RSpace::update making first copy of params", Logger::Trace);
      m_working_set.resize(parameters.size());
      std::iota(begin(m_working_set), end(m_working_set), size_t{0});
      ind_last_param_to_new = m_working_set;
    } else {
      m_logger->msg("RSpace::update updating params with new values ", Logger::Trace);
      ind_last_param_to_new = rspace::assign_last_parameters_to_new(m_last_params, parameters, m_handlers->qr());
    }
    m_params.clear();
    m_actions.clear();
    for (auto i : ind_last_param_to_new) {
      m_params.emplace_back(parameters.at(i));
      m_actions.emplace_back(actions.at(i));
    }
    for (size_t i = 0; i < m_params.size(); ++i) {
      auto norm = std::sqrt(m_handlers->rr().dot(m_params[i], m_params[i]));
      m_logger->msg("RSpace::update param index i = " + std::to_string(i) + ", 1./norm = " + Logger::scientific(norm),
                    Logger::Debug);
      // FIXME What happens if norm is very large or very small?
      m_handlers->rr().scal(1.0 / norm, m_params[i]);
      m_handlers->rr().scal(1.0 / norm, m_actions[i]);
    }
    data[EqnData::S] = util::overlap(m_params, m_handlers->rr());
    data[EqnData::H] = util::overlap(m_params, m_actions, m_handlers->rr());
    if (m_logger->data_dump) {
      m_logger->msg("S = " + as_string(data[EqnData::S]), Logger::Info);
      m_logger->msg("H = " + as_string(data[EqnData::H]), Logger::Info);
    }
    assert(m_working_set.size() == m_params.size());
  }

  size_t size() { return m_params.size(); }

  /*!
   * @brief Returns dummy containers that can be used as intermediates. Some parameters must first be assigned.
   * @param n number of dummy containers required
   */
  auto& dummy(size_t n) const {
    assert(!m_params.empty() && "must add parameters to the RSpace first");
    while (m_dummy.size() < n) {
      m_dummy.emplace_back(m_handlers->rr().copy(m_params.front()));
    }
    return m_dummy;
  }

  //! Updates working set of vectors. @param working_vector_ind indices of params that are still part of the working set
  void update_working_set(const std::vector<size_t>& working_vector_ind) {
    assert(working_vector_ind.size() <= m_params.size());
    auto n_copy = std::min(m_last_params.size(), working_vector_ind.size());
    auto n_tot = working_vector_ind.size();
    for (size_t i = 0; i < n_copy; ++i) {
      m_handlers->qr().copy(m_last_params[i], m_params.at(working_vector_ind[i]));
      m_handlers->qr().copy(m_last_actions[i], m_actions.at(working_vector_ind[i]));
    }
    for (size_t i = n_copy; i < n_tot; ++i) {
      m_last_params.emplace_back(m_handlers->qr().copy(m_params.at(working_vector_ind[i])));
      m_last_actions.emplace_back(m_handlers->qr().copy(m_actions.at(working_vector_ind[i])));
    }
    m_last_params.resize(n_tot);
    m_last_actions.resize(n_tot);
    auto new_working_set = std::vector<size_t>{};
    for (auto ind : working_vector_ind) {
      new_working_set.emplace_back(m_working_set.at(ind));
    }
    m_logger->msg("RSpace::update_working_set old working set = ", begin(m_working_set), end(m_working_set),
                  Logger::Debug);
    m_logger->msg("RSpace::update_working_set new working set = ", begin(new_working_set), end(new_working_set),
                  Logger::Debug);
    m_working_set = new_working_set;
  }

  //! Returns list of root indices for each working vector. Each element corresponds to element in params.
  auto& working_set() const { return m_working_set; }

  auto& params() const { return m_params; }
  auto& params() { return m_params; }
  auto& actions() const { return m_actions; }
  auto& actions() { return m_actions; }
  auto& last_params() { return m_last_params; }
  auto& last_params() const { return m_last_params; }
  auto& last_actions() { return m_last_actions; }
  auto& last_actions() const { return m_last_actions; }

protected:
  std::shared_ptr<ArrayHandlers<R, Q, P>> m_handlers;
  std::shared_ptr<Logger> m_logger;
  std::vector<size_t> m_working_set; //!< working set of roots. Maps references of current params to starting roots
  VecRefR m_params;                  //!< solutions at this iteration forming the RSpace, mapped to root indices
  VecRefR m_actions;                 //!< action vector corresponding to m_params
  mutable std::vector<R> m_dummy;    //!< A dummy R vector which can be used as an intermediate
  std::vector<Q> m_last_params;      //!< parameters from previous iteration, mapped to root indices
  std::vector<Q> m_last_actions;     //!< actions from previous iteration, mapped to root indices
};

//! RSpace for LinearEquations solver
template <class R, class Q, class P>
class RSpaceLEq : public RSpace<R, Q, P> {
public:
  using RSpace<R, Q, P>::data;

  explicit RSpaceLEq(std::shared_ptr<ArrayHandlers<R, Q, P>> handlers) : RSpace<R, Q, P>(std::move(handlers)) {
    data = null_data<EqnData::H, EqnData::S, EqnData::rhs>;
  }

  void update(const std::vector<R>& parameters, const std::vector<R>& action, LinearEquations<R, Q, P>& solver) {
    RSpace<R, Q, P>::update_rspace(data, parameters, action, solver);
    // now update RHS vector
  }
};

} // namespace subspace
} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_RSPACE_H
