#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_RSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_RSPACE_H
#include <memory>
#include <molpro/linalg/itsolv/ArrayHandlers.h>
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/subspace/SubspaceData.h>
#include <molpro/linalg/itsolv/subspace/util.h>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace subspace {

namespace rspace {
//! Assigns new parameters to previous based on maximum of their overlap.
template <class R, class Q>
std::map<size_t, size_t> assign_new_parameters_to_last(const std::vector<R>& new_params,
                                                       const std::map<size_t, Q>& last_params,
                                                       array::ArrayHandler<Q, R>& handler) {
  using util::overlap;
  using util::wrap;
  assert(new_params.size() >= last_params.size());
  auto ov = overlap(wrap(last_params), wrap(new_params), handler);
  auto assign_param_to_last = std::map<size_t, size_t>{};
  while (!ov.empty()) {
    auto imax = std::max_element(begin(ov.data()), end(ov.data()));
    auto ind = ov.to_coord(distance(ov.data().begin(), imax));
    assign_param_to_last[ind.second] = ind.first;
    ov.remove_row(ind.first);
  }
  assert(assign_param_to_last.size() == last_params.size());
  return assign_param_to_last;
}

} // namespace rspace

//! Space for the working set of vectors
template <class R, class Q, class P>
class RSpace {
public:
  using MapRefR = std::map<size_t, std::reference_wrapper<R>>;

  //! Matrix and overlap data mapped to the subspace
  SubspaceData data = null_data<EqnData::H, EqnData::S>();

  explicit RSpace(std::shared_ptr<ArrayHandlers<R, Q, P>> handlers) : m_handlers(std::move(handlers)) {}

  void update(const std::vector<R>& parameters, const std::vector<R>& actions, IterativeSolver<R, Q, P>& solver) {
    auto assign_param_to_last = std::map<size_t, size_t>{};
    if (m_last_params.empty())
      for (size_t i = 0; i < parameters; ++i)
        assign_param_to_last[i] = i;
    else
      assign_param_to_last = rspace::assign_new_parameters_to_last(parameters, m_last_params, m_handlers->qr());
    m_working_set.clear();
    for (const auto& param_to_last : assign_param_to_last) {
      m_params[param_to_last.second] = parameters[param_to_last.first];
      m_actions[param_to_last.second] = actions[param_to_last.first];
      m_working_set.emplace_back(param_to_last.second);
    }
    for (size_t i = 0; i < m_params.size(); ++i) {
      auto norm = std::sqrt(m_handlers->rr().dot(m_params[i], m_params[i]));
      // FIXME What happens if norm is very large or very small?
      m_handlers->rr.scal(1.0 / norm, m_params[i]);
      m_handlers->rr.scal(1.0 / norm, m_actions[i]);
    }
    data[EqnData::S] = util::overlap(wrap(m_params), m_handlers->rr());
    data[EqnData::H] = util::overlap(wrap(m_params), wrap(m_actions), m_handlers->rr());
  }

  size_t size() { return m_params.size(); }

  auto& dummy() const {
    assert(!m_params.empty() && "must add parameters to the RSpace first");
    //    if (m_dummy.empty()){
    //      m_dummy.emplace_back(std::make_shared<R>(m_handlers->rr().copy(m_params.front())));
    //      m_dummy.emplace_back(std::make_shared<R>(m_handlers->rr().copy(m_actions.front())));
    //    }
    return *m_dummy;
  }

  auto& working_set() const { return m_working_set; }

  auto& params() const { return m_params; }
  auto& actions() const { return m_actions; }

  auto& last_params() { return m_last_params; }
  auto& last_params() const { return m_last_params; }
  auto& last_actions() { return m_last_actions; }
  auto& last_actions() const { return m_last_actions; }

protected:
  std::shared_ptr<ArrayHandlers<R, Q, P>> m_handlers;
  std::vector<size_t> m_working_set; //!< working set of roots
  MapRefR m_params;                  //!< solutions at this iteration forming the RSpace, mapped to root indices
  MapRefR m_actions;                 //!< action vector corresponding to m_params
  mutable std::vector<std::shared_ptr<R>> m_dummy; //!< A dummy R vector which can be used as an intermediate
  std::map<size_t, Q> m_last_params;               //!< parameters from previous iteration, mapped to root inidices
  std::map<size_t, Q> m_last_actions;              //!< actions from previous iteration, mapped to root inidices
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
