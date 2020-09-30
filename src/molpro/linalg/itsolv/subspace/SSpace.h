#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SSPACE_H
#include <molpro/linalg/itsolv/ArrayHandlers.h>
#include <molpro/linalg/itsolv/Logger.h>
#include <molpro/linalg/itsolv/subspace/SubspaceData.h>

#include <vector>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace subspace {

//! Space storing best set of solutions
template <class Rt, class Qt, class Pt, typename ST>
class SSpace {
public:
  using R = Rt;
  using Q = Qt;
  using P = Pt;
  using scalar_type = ST;

  //! Matrix and overlap data mapped to the subspace
  SubspaceData data = null_data<EqnData::H, EqnData::S>();

  explicit SSpace(std::shared_ptr<ArrayHandlers<R, Q, P>> handlers, std::shared_ptr<Logger> logger)
      : m_handlers(std::move(handlers)), m_logger(std::move(logger)) {}

  void update(int root, const R& param, const R& action, scalar_type error) {}

  const auto& params() const { return m_params; };
  auto& params() { return m_params; };
  const auto& actions() const { return m_actions; };
  auto& actions() { return m_actions; };
  const auto& errors() const { return m_errors; };

protected:
  std::shared_ptr<ArrayHandlers<R, Q, P>> m_handlers;
  std::shared_ptr<Logger> m_logger;
  std::vector<Q> m_params;
  std::vector<Q> m_actions;
  std::vector<scalar_type> m_errors;
};

} // namespace subspace
} // namespace itsolv
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SSPACE_H
