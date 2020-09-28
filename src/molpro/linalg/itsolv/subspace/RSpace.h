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

  // FIXME are we still tracking roots? Parameters map directly to the previous solutions, but we map solutions to roots
  /*! Takes a list of parameters and corresponding action vectors
   * @brief Updates the RSpace. New parameters are stored as d vectors. RSpace consists of previous solutions.
   *
   *
   * @param parameters new preconditioned parameters. Each parameter corresponds to a working set.
   * @param actions action of matrix on parameters
   * @param solver overlying solver
   */
  void update(std::vector<R>& parameters, std::vector<R>& actions, IterativeSolver<R, Q, P>& solver) {
    // FIXME reordering is unnecessary now
    if (m_dparams.empty()) {
      m_logger->msg("RSpace::update making first copy of params", Logger::Trace);
      m_working_set.resize(parameters.size());
      std::iota(begin(m_working_set), end(m_working_set), size_t{0});
      for (size_t i = 0; i < m_working_set.size(); ++i) {
        m_dparams.emplace_back(m_handlers->rr().copy(parameters.at(i)));
        m_dactions.emplace_back(m_handlers->rr().copy(actions.at(i)));
      }
    } else {
      m_logger->msg("RSpace::update updating dparams with new values ", Logger::Trace);
      for (size_t i = 0; i < m_working_set.size(); ++i) {
        m_handlers->rr().copy(m_dparams.at(i), parameters.at(i));
        m_handlers->rr().copy(m_dactions.at(i), actions.at(i));
      }
      m_dparams.resize(m_working_set.size());
      m_dactions.resize(m_working_set.size());
    }
    m_working_params.clear();
    m_working_actions.clear();
    for (size_t i = 0; i < m_working_set.size(); ++i) {
      m_working_params.emplace_back(parameters.at(i));
      m_working_actions.emplace_back(actions.at(i));
    }
    assert(m_working_set.size() == m_dparams.size());
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

  /*!
   * @brief Updates working set of vectors. RSpace is populated with current solutions
   * @param working_vector_ind
   */
  void update_working_set(const std::vector<size_t>& working_vector_ind) {
    assert(working_vector_ind.size() <= m_working_params.size());
    auto n_copy = std::min(m_params.size(), working_vector_ind.size());
    auto n_tot = working_vector_ind.size();
    for (size_t i = 0; i < n_copy; ++i) {
      m_handlers->rr().copy(m_params[i], m_working_params.at(working_vector_ind[i]));
      m_handlers->rr().copy(m_actions[i], m_working_actions.at(working_vector_ind[i]));
    }
    for (size_t i = n_copy; i < n_tot; ++i) {
      m_params.emplace_back(m_handlers->rr().copy(m_working_params.at(working_vector_ind[i])));
      m_actions.emplace_back(m_handlers->rr().copy(m_working_actions.at(working_vector_ind[i])));
    }
    m_params.resize(n_tot);
    m_actions.resize(n_tot);
    auto new_working_set = std::vector<size_t>{};
    for (auto ind : working_vector_ind) {
      new_working_set.emplace_back(m_working_set.at(ind));
    }
    m_logger->msg("RSpace::update_working_set old working set = ", begin(m_working_set), end(m_working_set),
                  Logger::Debug);
    m_logger->msg("RSpace::update_working_set new working set = ", begin(new_working_set), end(new_working_set),
                  Logger::Debug);
    m_working_set = new_working_set;
    data[EqnData::S] = util::overlap(util::wrap(m_params), m_handlers->rr());
    data[EqnData::H] = util::overlap(util::wrap(m_params), util::wrap(m_actions), m_handlers->rr());
    if (m_logger->data_dump) {
      m_logger->msg("S = " + as_string(data[EqnData::S]), Logger::Info);
      m_logger->msg("H = " + as_string(data[EqnData::H]), Logger::Info);
    }
  }

  //! Returns list of root indices for each working vector. Each element corresponds to element in params.
  auto& working_set() const { return m_working_set; }

  auto& params() const { return m_params; }
  auto& params() { return m_params; }
  auto& actions() const { return m_actions; }
  auto& actions() { return m_actions; }
  auto& dparams() { return m_dparams; }
  auto& dparams() const { return m_dparams; }
  auto& dactions() { return m_dactions; }
  auto& dactions() const { return m_dactions; }

protected:
  std::shared_ptr<ArrayHandlers<R, Q, P>> m_handlers;
  std::shared_ptr<Logger> m_logger;
  std::vector<size_t> m_working_set; //!< working set of roots. Maps references of current params to starting roots
  std::vector<R> m_params;           //!< solutions at this iteration forming the RSpace, mapped to root indices
  std::vector<R> m_actions;          //!< action vector corresponding to m_params
  VecRefR m_working_params;          //!< references to input parameters
  VecRefR m_working_actions;         //!< references to input actions
  mutable std::vector<R> m_dummy;    //!< A dummy R vector which can be used as an intermediate
  std::vector<Q> m_dparams;          //!< parameters from previous iteration, mapped to root indices
  std::vector<Q> m_dactions;         //!< actions from previous iteration, mapped to root indices
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
