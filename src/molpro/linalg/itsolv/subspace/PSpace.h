#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_PSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_PSPACE_H
#include <molpro/linalg/array/ArrayHandler.h>
#include <molpro/linalg/itsolv/subspace/SubspaceData.h>
#include <molpro/linalg/itsolv/wrap.h>

namespace molpro::linalg::itsolv::subspace {

template <class Rt, class Pt>
class PSpace {
public:
  using R = Rt;
  using P = Pt;

  void update(const CVecRef<P>& params, array::ArrayHandler<P, P>& handler) {
    for (const auto& p : params)
      m_params.emplace_back(handler.copy(p));
  }

  CVecRef<P> params() const { return cwrap(m_params); }
  CVecRef<P> cparams() const { return params(); }
  VecRef<P> params() { return wrap(m_params); }
  CVecRef<P> actions() const { return {}; }
  CVecRef<P> cactions() const { return actions(); }
  VecRef<P> actions() { return {}; }

  size_t size() const { return m_params.size(); }

  void erase(size_t i) { m_params.erase(m_params.begin() + i); }

private:
  std::vector<P> m_params;
};

} // namespace molpro::linalg::itsolv::subspace

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_PSPACE_H
