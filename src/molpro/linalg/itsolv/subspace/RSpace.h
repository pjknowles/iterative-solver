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

  /*!
   * @brief Updates the RSpace with the new parameters.
   *
   *
   * @param parameters new preconditioned parameters. Each parameter corresponds to a working set.
   * @param actions action of matrix on parameters
   * @param solver overlying solver
   */
  void update(std::vector<R>& parameters, std::vector<R>& actions, IterativeSolver<R, Q, P>& solver) {
    m_logger->msg("RSpace::update", Logger::Trace);
    assert(parameters.size() == actions.size());
    clear();
    std::copy_n(begin(parameters), solver.working_set().size(), std::back_inserter(m_params));
    std::copy_n(begin(actions), solver.working_set().size(), std::back_inserter(m_actions));
    data[EqnData::S] = util::overlap(m_params, m_handlers->rr());
    data[EqnData::H] = util::overlap(m_params, m_actions, m_handlers->rr());
    if (m_logger->data_dump) {
      m_logger->msg("Srr = " + as_string(data[EqnData::S]), Logger::Info);
      m_logger->msg("Hrr = " + as_string(data[EqnData::H]), Logger::Info);
    }
  }

  //! Number of parameters in the subspace
  size_t size() const { return m_params.size(); }

  //! Remove all parameters from subspace
  void clear() {
    m_params.clear();
    m_actions.clear();
    for (auto& datum : data) {
      datum.second.clear();
    }
  }

  //! Remove parameter at index i
  void erase(size_t i) {
    m_params.erase(std::next(m_params.begin(), i));
    m_actions.erase(std::next(m_actions.begin(), i));
    for (auto& datum : data) {
      datum.second.remove_row_col(i, i);
    }
  }

  /*!
   * @brief Returns dummy containers that can be used as intermediates. Some parameters must first be assigned.
   * @param n number of dummy containers required
   */
  auto& dummy(size_t n) const {
    assert((!m_params.empty()) && "must add parameters to the RSpace first");
    while (m_dummy.size() < n) {
      m_dummy.emplace_back(m_handlers->rr().copy(m_params.front()));
    }
    return m_dummy;
  }

  auto& params() const { return m_params; }
  auto& params() { return m_params; }
  auto& actions() const { return m_actions; }
  auto& actions() { return m_actions; }

protected:
  std::shared_ptr<ArrayHandlers<R, Q, P>> m_handlers;
  std::shared_ptr<Logger> m_logger;
  VecRefR m_params;               //!< parameters from user
  VecRefR m_actions;              //!< actions corresponding to m_params
  mutable std::vector<R> m_dummy; //!< A dummy R vector which can be used as an intermediate
};

} // namespace subspace
} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_RSPACE_H
