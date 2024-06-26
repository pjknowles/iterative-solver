#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_IXSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_IXSPACE_H
#include <molpro/linalg/array/ArrayHandler.h>
#include <molpro/linalg/array/Span.h>
#include <molpro/linalg/itsolv/subspace/Dimensions.h>
#include <molpro/linalg/itsolv/subspace/SubspaceData.h>
#include <molpro/linalg/itsolv/wrap.h>

namespace molpro::linalg::itsolv::subspace {

//! Full subspace
template <class RT, class QT, class PT>
class IXSpace {
public:
  using R = RT;
  using Q = QT;
  using P = PT;
  using value_type = typename array::ArrayHandler<R, R>::value_type;
  using value_type_abs = typename array::ArrayHandler<R, R>::value_type_abs;
  IXSpace() = default;
  virtual ~IXSpace() = default;
  SubspaceData data; //!< Equation data in the subspace

  //! Number of vectors forming the subspace
  size_t size() const { return dimensions().nX; }

  //! Removes parameter i from the full subspace
  virtual void erase(size_t i) = 0;
  //! Removes parameter i from Q subspace
  virtual void eraseq(size_t i) = 0;
  //! Removes parameter i from P subspace
  virtual void erasep(size_t i) = 0;
  //! Removes parameter i from D subspace
  virtual void erased(size_t i) = 0;

  //! Adds parameters to the P space
  virtual void update_pspace(const CVecRef<P>& params, const array::Span<value_type>& pp_action_matrix) = 0;

  //! Adds parameters to the Q space
  virtual void update_qspace(const CVecRef<R>& params, const CVecRef<R>& actions) = 0;

  //! Updates D space with the new parameters
  virtual void update_dspace(VecRef<Q>& params, VecRef<Q>& actions) = 0;

  virtual VecRef<P> paramsp() = 0;
  virtual VecRef<Q> paramsq() = 0;
  virtual VecRef<Q> actionsq() = 0;
  virtual VecRef<Q> paramsd() = 0;
  virtual VecRef<Q> actionsd() = 0;

  virtual CVecRef<P> paramsp() const = 0;
  virtual CVecRef<Q> paramsq() const = 0;
  virtual CVecRef<Q> actionsq() const = 0;
  virtual CVecRef<Q> paramsd() const = 0;
  virtual CVecRef<Q> actionsd() const = 0;

  virtual CVecRef<P> cparamsp() const = 0;
  virtual CVecRef<Q> cparamsq() const = 0;
  virtual CVecRef<Q> cactionsq() const = 0;
  virtual CVecRef<Q> cparamsd() const = 0;
  virtual CVecRef<Q> cactionsd() const = 0;

  virtual const Dimensions& dimensions() const = 0;
};

} // namespace molpro::linalg::itsolv::subspace

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_IXSPACE_H