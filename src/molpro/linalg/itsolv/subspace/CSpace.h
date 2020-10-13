#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_CSPACE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_CSPACE_H
#include <molpro/linalg/itsolv/ArrayHandlers.h>
#include <molpro/linalg/itsolv/Logger.h>
#include <molpro/linalg/itsolv/subspace/Dimensions.h>
#include <molpro/linalg/itsolv/subspace/SubspaceData.h>
#include <molpro/linalg/itsolv/wrap.h>

#include <vector>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace subspace {
/*!
 * @brief Space storing the complement of Q necessary to reconstruct previous solutions.
 *
 * It consists of linear combinations of deleted Q parameters and maximum size of the space is the number of solutions.
 * As such, C space accumulates round off error and eventually actions will not map to corresponding parameters
 * leading to errors in the subspace. Thus C space should be periodically reset.
 *
 * Resetting C space
 * -----------------
 *  Use C space as new parameters, calculate their exact action and add to Q.
 *
 *  Cons: size of C space can be greater than R, so what will happen to the rest of C space?
 *
 * This can be done over multiple iterations, as long as no more Q parameters are removed. Thus firstly remove nC
 * parameters from Q, than start moving C parameters to R params in propose_rspace, until C space is empty. At that
 * point all of C space was moved into Q and there were no modifications, since size of Q space was within limit and
 * residual was not used. It should be possible to have a mixed update where new parameters are orthogonalised residuals
 * and C params.
 */
template <class Rt, class Qt, class Pt, typename ST>
class CSpace {
public:
  using R = Rt;
  using Q = Qt;
  using P = Pt;
  using scalar_type = ST;

  explicit CSpace(std::shared_ptr<ArrayHandlers<R, Q, P>> handlers, std::shared_ptr<Logger> logger)
      : m_handlers(std::move(handlers)), m_logger(std::move(logger)) {}

  //! Adds a solution to the space, overwriting if a solution for the root already exists, and construct equation data
  void update(const unsigned int root, const R& param, const R& action, scalar_type error) {
    auto it = m_params.find(root);
    if (it == m_params.end()) {
      m_params.emplace(root, m_handlers->qr().copy(param));
      m_actions.emplace(root, m_handlers->qr().copy(action));
    } else {
      m_handlers->qr().copy(m_params.at(root), param);
      m_handlers->qr().copy(m_actions.at(root), action);
    }
  }

  //! Adds solutions to the space, overwriting if a solution for given root already exists, and construct equation data
  void update(const std::vector<unsigned int>& roots, const CVecRef<R>& params, const CVecRef<R>& actions,
              const std::vector<scalar_type>& errors) {
    for (size_t i = 0; i < roots.size(); ++i) {
      update(roots.at(i), params.at(i), actions.at(i), errors.at(i));
    }
  }

  size_t size() const { return m_params.size(); }

  //! Removes all vectors and data
  void clear() {
    m_params.clear();
    m_actions.clear();
  }

  //! Erases parameter i. @param i index in the current space
  void erase(size_t i) {
    assert(m_params.size() > i);
    auto erase_at_i = [i](auto& v) { v.erase(std::next(begin(v), i)); };
    erase_at_i(m_params);
    erase_at_i(m_actions);
  }

  CVecRef<Q> params() const { return wrap(m_params); };
  CVecRef<Q> cparams() const { return params(); };
  VecRef<Q> params() { return wrap(m_params); };
  CVecRef<Q> actions() const { return wrap(m_actions); };
  CVecRef<Q> cactions() const { return actions(); };
  VecRef<Q> actions() { return wrap(m_actions); };

protected:
  std::shared_ptr<ArrayHandlers<R, Q, P>> m_handlers;
  std::shared_ptr<Logger> m_logger;
  std::map<unsigned int, Q> m_params;  //! parameters for root index
  std::map<unsigned int, Q> m_actions; //! actions for root index
};

} // namespace subspace
} // namespace itsolv
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_CSPACE_H
