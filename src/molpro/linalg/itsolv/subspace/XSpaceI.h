#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_XSPACEI_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_XSPACEI_H
#include <molpro/linalg/itsolv/subspace/Dimensions.h>
#include <molpro/linalg/itsolv/subspace/SubspaceData.h>
#include <molpro/linalg/itsolv/wrap.h>

namespace molpro::linalg::itsolv::subspace {

//! Full subspace
template <class RT, class QT, class PT>
class XSpaceI {
public:
  using R = RT;
  using Q = QT;
  using P = PT;
  using value_type = typename array::ArrayHandler<R, R>::value_type;
  using value_type_abs = typename array::ArrayHandler<R, R>::value_type_abs;
  XSpaceI() = default;
  virtual ~XSpaceI() = default;
  SubspaceData data; //!< Equation data in the subspace

  //! Number of vectors forming the subspace
  size_t size() { return dimensions().nX; }

  //! Removes parameter i from the full subspace
  virtual void erase(size_t i) = 0;
  //! Removes parameter i from Q subspace
  virtual void eraseq(size_t i) = 0;
  //! Removes parameter i from P subspace
  virtual void erasep(size_t i) = 0;
  //! Removes parameter i from D subspace
  virtual void erased(size_t i) = 0;

  //! Adds parameters to the Q space
  virtual void update_pspace() {}

  //! Adds parameters to the Q space
  virtual void update_qspace(const CVecRef<R>& params, const CVecRef<R>& actions) = 0;

  //! Updates D space with the new parameters
  virtual void update_dspace(VecRef<R>& params, VecRef<R>& actions, const Matrix<value_type>& lin_trans_only_R) = 0;

  //! Completes construction of D space action by adding contribution from the R space, and constructs equation data
  //! blocks corresponding to D space
  virtual void complete_dspace_action(const CVecRef<R>& actions) = 0;

  virtual VecRef<P> paramsp() = 0;
  virtual VecRef<P> actionsp() = 0;
  virtual VecRef<Q> paramsq() = 0;
  virtual VecRef<Q> actionsq() = 0;
  virtual VecRef<Q> paramsd() = 0;
  virtual VecRef<Q> actionsd() = 0;

  virtual CVecRef<P> paramsp() const = 0;
  virtual CVecRef<P> actionsp() const = 0;
  virtual CVecRef<Q> paramsq() const = 0;
  virtual CVecRef<Q> actionsq() const = 0;
  virtual CVecRef<Q> paramsd() const = 0;
  virtual CVecRef<Q> actionsd() const = 0;

  virtual CVecRef<P> cparamsp() const = 0;
  virtual CVecRef<P> cactionsp() const = 0;
  virtual CVecRef<Q> cparamsq() const = 0;
  virtual CVecRef<Q> cactionsq() const = 0;
  virtual CVecRef<Q> cparamsd() const = 0;
  virtual CVecRef<Q> cactionsd() const = 0;

  virtual const xspace::Dimensions& dimensions() const = 0;
};

} // namespace molpro::linalg::itsolv::subspace

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_XSPACEI_H