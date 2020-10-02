#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_CSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_CSPACE_H
#include <molpro/linalg/itsolv/ArrayHandlers.h>
#include <molpro/linalg/itsolv/Logger.h>
#include <molpro/linalg/itsolv/subspace/SubspaceData.h>
#include <molpro/linalg/itsolv/wrap.h>

#include <vector>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace subspace {

//! Space storing the best set of solutions and their errors
template <class Rt, class Qt, class Pt, typename ST>
class CSpace {
public:
  using R = Rt;
  using Q = Qt;
  using P = Pt;
  using scalar_type = ST;

  //! Matrix and overlap data mapped to the subspace
  SubspaceData data = null_data<EqnData::H, EqnData::S>();

  explicit CSpace(std::shared_ptr<ArrayHandlers<R, Q, P>> handlers, std::shared_ptr<Logger> logger)
      : m_handlers(std::move(handlers)), m_logger(std::move(logger)) {}

  //! Adds a solution to the space, overwriting if a solution for the root already exists
  void update(unsigned int root, const R& param, const R& action, scalar_type error) {
    if (root >= m_params.size()) {
      m_params.emplace_back(m_handlers->qr().copy(param));
      m_params.emplace_back(m_handlers->qr().copy(action));
      m_errors.emplace_back(error);
    } else {
      m_handlers->qr().copy(m_params.at(root), param);
      m_handlers->qr().copy(m_params.at(root), action);
      m_errors.at(root) = error;
    }
  }

  //! Adds solutions to the space, overwriting if a solution for given root already exists
  void update(const std::vector<unsigned int>& roots, const std::vector<R>& params, const std::vector<R>& actions,
              const std::vector<scalar_type>& errors) {
    for (size_t i = 0; i < roots.size(); ++i) {
      update(roots.at(i), params.at(i), actions.at(i), errors.at(i));
    }
  }

  void set_error(unsigned int root, scalar_type error) {}
  void set_error(const std::vector<unsigned int>& roots, const std::vector<scalar_type>& errors) {}

  size_t size() const { return m_params.size(); }

  //! Removes all vectors and data
  void clear() {
    m_params.clear();
    m_actions.clear();
    m_errors.clear();
    for (auto d : data) {
      d.second.clear();
    }
  }

  //! Erases parameter i. @param i index in the current space
  void erase(size_t i) {
    assert(m_params.size() > i);
    auto it = std::next(begin(m_params), i);
    auto erase_at_i = [i](auto& v) { v.erase(std::next(begin(v), i)); };
    for (auto& v : {m_params, m_actions, m_errors})
      erase_at_i(v);
    for (auto d : {EqnData::H, EqnData::S}) {
      data.at(d).remove_row_col(i, i);
    }
  }

  CVecRef<Q> params() const { return wrap(m_params); };
  CVecRef<Q> cparams() const { return params(); };
  VecRef<Q> params() { return wrap(m_params); };
  CVecRef<Q> actions() const { return wrap(m_actions); };
  CVecRef<Q> cactions() const { return actions(); };
  VecRef<Q> actions() { return wrap(m_actions); };
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
