#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_CSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_CSPACE_H
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
class CSpace {
public:
  using R = Rt;
  using Q = Qt;
  using P = Pt;
  using VecRefQ = std::vector<std::reference_wrapper<Q>>;
  using CVecRefQ = std::vector<std::reference_wrapper<const Q>>;
  using scalar_type = ST;

  //! Matrix and overlap data mapped to the subspace
  SubspaceData data = null_data<EqnData::H, EqnData::S>();

  explicit CSpace(std::shared_ptr<ArrayHandlers<R, Q, P>> handlers, std::shared_ptr<Logger> logger)
      : m_handlers(std::move(handlers)), m_logger(std::move(logger)) {}

  void update(size_t root, const R& param, const R& action, scalar_type error) {}
  void update(const std::vector<size_t>& roots, const std::vector<R>& params, const std::vector<R>& actions,
              const std::vector<scalar_type>& errors) {}

  void set_error(int root, scalar_type error) {}
  void set_error(const std::vector<size_t>& roots, const std::vector<scalar_type>& errors) {}

  size_t size() const { return m_params.size(); }

  CVecRefQ params() const {
    auto params = CVecRefQ{};
    std::transform(begin(m_params), end(m_params), std::back_inserter(params), [](auto& p) { return std::ref(p); });
    return params;
  };

  VecRefQ params() {
    auto params = VecRefQ{};
    std::transform(begin(m_params), end(m_params), std::back_inserter(params), [](auto& p) { return std::ref(p); });
    return params;
  };

  CVecRefQ actions() const {
    auto actions = CVecRefQ{};
    std::transform(begin(m_actions), end(m_actions), std::back_inserter(actions), [](auto& p) { return std::ref(p); });
    return actions;
  };

  VecRefQ actions() {
    auto actions = VecRefQ{};
    std::transform(begin(m_actions), end(m_actions), std::back_inserter(actions), [](auto& p) { return std::ref(p); });
    return actions;
  };

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
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_CSPACE_H
