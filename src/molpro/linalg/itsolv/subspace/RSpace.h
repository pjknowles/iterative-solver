#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_RSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_RSPACE_H
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/subspace/SubspaceData.h>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace detail {

namespace rspace {

template <class R, class Q, class P>
void update_rspace(SubspaceData& rs, const std::vector<R>& parameters, const std::vector<R>& action);

template <class R, class Q, class P>
void update_rspace_rhs(SubspaceData& rs, const std::vector<R>& parameters, const std::vector<R>& action,
                       LinearEigensystem<R, Q, P>& solver);
} // namespace rspace

//!
template <class R, class Q>
class RSpace {
public:
  using MapRefR = std::map<size_t, std::reference_wrapper<R>>;

  SubspaceData subspace = null_data<EqnData::H, EqnData::S>();

  template <class P>
  void update(const std::vector<R>& parameters, const std::vector<R>& action, IterativeSolver<R, Q, P>& solver) {
    //    detail::update_rspace(subspace, parameters, action);
    auto working_set = solver.working_set();
    // normalise parameters
    // overlap matrix with previous parameters
    // get parameter order from the overlap
    // save references to parameters in that order
    // construct subspace
  }

  size_t size() { return subspace.at(EqnData::H).rows(); }

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
  std::vector<size_t> m_working_set; //!< working set of roots
  MapRefR m_params;                  //!< solutions at this iteration forming the RSpace, mapped to root indices
  MapRefR m_actions;                 //!< action vector corresponding to m_params
  mutable std::vector<std::shared_ptr<R>> m_dummy; //!< A dummy R vector which can be used as an intermediate
  std::map<size_t, Q> m_last_params;               //!< parameters from previous iteration, mapped to root inidices
  std::map<size_t, Q> m_last_actions;              //!< actions from previous iteration, mapped to root inidices
};

//! RSpace for LinearEquations solver
template <class R, class Q>
class RSpaceLEq : public RSpace<R, Q> {
public:
  using RSpace<R, Q>::subspace;
  RSpaceLEq() : RSpace<R, Q>() { subspace = null_data<EqnData::H, EqnData::S, EqnData::rhs>; }

  template <class P>
  void update(const std::vector<R>& parameters, const std::vector<R>& action, LinearEquations<R, Q, P>& solver) {
    RSpace<R, Q>::update_rspace(subspace, parameters, action, solver);
    // now update RHS vector
  }
};

} // namespace detail
} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_RSPACE_H
