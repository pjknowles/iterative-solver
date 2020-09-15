#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_RSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_RSPACE_H
#include <functional>
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
//! Assigns new parameters to previous based on maximum overlap.
template <class R, class Q>
std::vector<size_t> assign_new_parameters_to_last(const std::vector<R>& new_params, const std::vector<Q>& last_params,
                                                  array::ArrayHandler<Q, R>& handler) {
  using util::overlap;
  using util::wrap;
  assert(new_params.size() >= last_params.size());
  auto ov = overlap(wrap(last_params), wrap(new_params), handler);
  auto assign_param_to_last = std::vector<size_t>(last_params.size());
  for (size_t i = 0; i < ov.rows(); ++i) {
    auto imax = std::max_element(begin(ov.data()), end(ov.data()));
    auto ind = ov.to_coord(distance(ov.data().begin(), imax));
    assign_param_to_last[ind.first] = ind.second;
    ov.remove_row(ind.first);
  }
  return assign_param_to_last;
}

} // namespace rspace

//! Space for the working set of vectors
template <class R, class Q, class P>
class RSpace {
public:
  using VecRefR = std::vector<std::reference_wrapper<R>>;

  //! Matrix and overlap data mapped to the subspace
  SubspaceData data = null_data<EqnData::H, EqnData::S>();

  explicit RSpace(std::shared_ptr<ArrayHandlers<R, Q, P>> handlers) : m_handlers(std::move(handlers)) {}

  void update(std::vector<R>& parameters, std::vector<R>& actions, IterativeSolver<R, Q, P>& solver) {
    auto assign_param_to_last = std::vector<size_t>{};
    if (m_last_params.empty()) {
      m_working_set.resize(parameters.size());
      std::iota(begin(m_working_set), end(m_working_set), size_t{0});
      assign_param_to_last = m_working_set;
    } else
      assign_param_to_last = rspace::assign_new_parameters_to_last(parameters, m_last_params, m_handlers->qr());
    m_params.clear();
    m_actions.clear();
    for (const auto& param_to_last : assign_param_to_last) {
      m_params.emplace_back(parameters[param_to_last]);
      m_actions.emplace_back(actions[param_to_last]);
    }
    for (size_t i = 0; i < m_params.size(); ++i) {
      auto norm = std::sqrt(m_handlers->rr().dot(m_params[i], m_params[i]));
      // FIXME What happens if norm is very large or very small?
      m_handlers->rr().scal(1.0 / norm, m_params[i]);
      m_handlers->rr().scal(1.0 / norm, m_actions[i]);
    }
    data[EqnData::S] = util::overlap(m_params, m_handlers->rr());
    data[EqnData::H] = util::overlap(m_params, m_actions, m_handlers->rr());
    assert(m_working_set.size() == m_params.size());
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

  //! Returns list of working set of roots. Each element corresponds to element in params.
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
  std::vector<size_t> m_working_set; //!< working set of roots. Maps references of current params to starting roots
  VecRefR m_params;                  //!< solutions at this iteration forming the RSpace, mapped to root indices
  VecRefR m_actions;                 //!< action vector corresponding to m_params
  mutable std::vector<std::shared_ptr<R>> m_dummy; //!< A dummy R vector which can be used as an intermediate
  std::vector<Q> m_last_params;                    //!< parameters from previous iteration, mapped to root indices
  std::vector<Q> m_last_actions;                   //!< actions from previous iteration, mapped to root indices
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
