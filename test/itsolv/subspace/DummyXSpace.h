#ifndef LINEARALGEBRA_TEST_ITSOLV_SUBSPACE_DUMMYXSPACE_H
#define LINEARALGEBRA_TEST_ITSOLV_SUBSPACE_DUMMYXSPACE_H

#include <molpro/linalg/itsolv/subspace/Dimensions.h>
#include <molpro/linalg/itsolv/subspace/XSpaceI.h>

using molpro::linalg::itsolv::CVecRef;
using molpro::linalg::itsolv::VecRef;
using molpro::linalg::itsolv::subspace::Matrix;

namespace {
template <class R, class Q, class P>
struct DummyXSpace : public molpro::linalg::itsolv::subspace::XSpaceI<R, Q, P> {
  using typename molpro::linalg::itsolv::subspace::XSpaceI<R, Q, P>::value_type;
  using typename molpro::linalg::itsolv::subspace::XSpaceI<R, Q, P>::value_type_abs;
  //! Removes parameter i from the full subspace
  void erase(size_t i) override {}
  //! Removes parameter i from Q subspace
  void eraseq(size_t i) override {}
  //! Removes parameter i from P subspace
  void erasep(size_t i) override {}
  //! Removes parameter i from D subspace
  void erased(size_t i) override {}

  //! Adds parameters to the Q space
  void update_pspace(const CVecRef<P>& params,
                     const molpro::linalg::array::Span<value_type>& pp_action_matrix) override {}

  //! Adds parameters to the Q space
  void update_qspace(const CVecRef<R>& params, const CVecRef<R>& actions) override {}

  //! Updates D space with the new parameters
  void update_dspace(VecRef<Q>& params, VecRef<Q>& actions) override {}

  VecRef<P> paramsp() override { return VecRef<P>{}; }
  VecRef<Q> paramsq() override { return VecRef<Q>{}; }
  VecRef<Q> actionsq() override { return VecRef<Q>{}; }
  VecRef<Q> paramsd() override { return VecRef<Q>{}; }
  VecRef<Q> actionsd() override { return VecRef<Q>{}; }

  CVecRef<P> paramsp() const override { return CVecRef<P>{}; }
  CVecRef<Q> paramsq() const override { return CVecRef<Q>{}; }
  CVecRef<Q> actionsq() const override { return CVecRef<Q>{}; }
  CVecRef<Q> paramsd() const override { return CVecRef<Q>{}; }
  CVecRef<Q> actionsd() const override { return CVecRef<Q>{}; }

  CVecRef<P> cparamsp() const override { return CVecRef<P>{}; }
  CVecRef<Q> cparamsq() const override { return CVecRef<Q>{}; }
  CVecRef<Q> cactionsq() const override { return CVecRef<Q>{}; }
  CVecRef<Q> cparamsd() const override { return CVecRef<Q>{}; }
  CVecRef<Q> cactionsd() const override { return CVecRef<Q>{}; }

  const molpro::linalg::itsolv::subspace::Dimensions& dimensions() const override { return dims; }

  molpro::linalg::itsolv::subspace::Dimensions dims;
};
} // namespace

#endif // LINEARALGEBRA_TEST_ITSOLV_SUBSPACE_DUMMYXSPACE_H
