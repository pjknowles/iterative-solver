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
  //! Removes parameter i from C subspace
  virtual void erasec(size_t i) = 0;

  //! Adds parameters to the Q space
  virtual void update_pspace() {}

  //! Adds parameters to the Q space
  virtual void update_qspace(const CVecRef<R>& params, const CVecRef<R>& actions) = 0;

  //! Adds solutions to the C space
  virtual void update_cspace(const std::vector<unsigned int>& roots, const CVecRef<R>& params,
                             const CVecRef<R>& actions, const std::vector<value_type>& errors) = 0;

  virtual VecRef<P> paramsp() = 0;
  virtual VecRef<P> actionsp() = 0;
  virtual VecRef<Q> paramsq() = 0;
  virtual VecRef<Q> actionsq() = 0;
  virtual VecRef<Q> paramsc() = 0;
  virtual VecRef<Q> actionsc() = 0;

  virtual CVecRef<P> paramsp() const = 0;
  virtual CVecRef<P> actionsp() const = 0;
  virtual CVecRef<Q> paramsq() const = 0;
  virtual CVecRef<Q> actionsq() const = 0;
  virtual CVecRef<Q> paramsc() const = 0;
  virtual CVecRef<Q> actionsc() const = 0;

  virtual CVecRef<P> cparamsp() const = 0;
  virtual CVecRef<P> cactionsp() const = 0;
  virtual CVecRef<Q> cparamsq() const = 0;
  virtual CVecRef<Q> cactionsq() const = 0;
  virtual CVecRef<Q> cparamsc() const = 0;
  virtual CVecRef<Q> cactionsc() const = 0;

  virtual const xspace::Dimensions& dimensions() const = 0;
};

} // namespace molpro::linalg::itsolv::subspace

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_XSPACEI_H