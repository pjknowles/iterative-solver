#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_Q_H_
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_Q_H_
#include <cmath>
#include <map>
#include <memory>
#include <molpro/linalg/iterativesolver/ArrayHandlers.h>
#include <molpro/linalg/iterativesolver/P.h>
//#include <molpro/iostream.h>
#include <vector>

namespace molpro {
namespace linalg {
namespace iterativesolver {

/*!
 * @brief A class holding a Q space
 * @tparam Rvector In-memory vectors. New instances will not be created.
 * @tparam Qvector Out-of memory vectors to be used for storing the space.
 * @tparam Pvector Containers that specify a P-space member
 */
template <class Rvector, class Qvector, class Pvector, typename scalar_type>
class Q {
  std::shared_ptr<ArrayHandlers<Rvector, Qvector, Pvector>> m_handlers;
  using key_type = int;
  bool m_hermitian;
  std::map<key_type, std::map<key_type, scalar_type>> m_metric;
  std::map<key_type, std::map<key_type, scalar_type>> m_action;
  std::map<key_type, std::map<key_type, scalar_type>> m_action_action;
  std::map<key_type, std::vector<scalar_type>> m_rhs;
  key_type m_index = 0;
  std::map<key_type, Qvector> m_vectors;
  std::map<key_type, Qvector> m_actions;
  const P<Pvector>& m_pspace;
  std::map<key_type, std::vector<scalar_type>> m_metric_pspace;
  std::map<key_type, std::vector<scalar_type>> m_action_pspace;
  std::vector<key_type> m_keys;
  std::map<key_type, typename Rvector::value_type> m_scale_factors, m_diff_factors;

public:
  Q(const P<Pvector>& pspace, bool hermitian, const std::shared_ptr<ArrayHandlers<Rvector, Qvector, Pvector>>& handlers)
      : m_hermitian(hermitian), m_pspace(pspace), m_handlers(handlers) {}

  const scalar_type& metric(size_t i, size_t j) const { return m_metric.at(m_keys[i]).at(m_keys[j]); }

  const scalar_type& action(size_t i, size_t j) const { return m_action.at(m_keys[i]).at(m_keys[j]); }
  const scalar_type& action_action(size_t i, size_t j) const { return m_action_action.at(m_keys[i]).at(m_keys[j]); }

  const std::vector<scalar_type>& rhs(size_t i) const { return m_rhs.at(m_keys[i]); }

  const std::vector<scalar_type>& metric_pspace(size_t i) const { return m_metric_pspace.at(m_keys[i]); }
  const std::vector<scalar_type>& action_pspace(size_t i) const { return m_action_pspace.at(m_keys[i]); }

  size_t size() const { return m_vectors.size(); }

  const Qvector& operator[](size_t i) const { return m_vectors.at(m_keys[i]); }
  const Qvector& action(size_t i) const { return m_actions.at(m_keys[i]); }

  typename Rvector::value_type scale_factor(size_t i) const { return m_scale_factors.at(m_keys[i]); }

public:
  /*!
   * @brief Obtain all of the keys that index vectors in the Q space.
   * If the Q space is modified such that vectors change their position, the key of the vector will never change.
   * @return
   */
  std::vector<key_type> keys() const {
    std::vector<key_type> result;
    for (const auto& vi : m_vectors)
      result.push_back(vi.first);
    return result;
  }

public:
  /*!
   * @brief  Assert or test whether the underlying kernel matrix in linear problems is hermitian.
   * @param hermitian The new state
   * @return The previous state
   */
  bool hermitian(bool hermitian = true) {
    auto old = m_hermitian;
    m_hermitian = hermitian;
    return old;
  }

  /*!
   * @brief Add a new vector to the Q space. Also compute and store the new elements of the QQ overlap and action
   * matrices, and overlap and interaction with P space.
   * @param vector
   * @param action
   * @param rhs
   * @param resres If true, action matrix will be action.action instead of vector.action
   */
  void add(const Rvector& vector, const Rvector& action, const std::vector<Qvector>& rhs, bool resres = false) {
    for (const auto& vi : (resres ? m_actions : m_vectors)) {
      m_metric[m_index][vi.first] = m_metric[vi.first][m_index] = m_handlers->rq().dot(vector, vi.second);
      m_action[vi.first][m_index] = m_handlers->rq().dot(action, vi.second);
    }
    for (const auto& vi : m_actions) {
      m_action_action[m_index][vi.first] = m_action_action[vi.first][m_index] = m_handlers->rq().dot(action, vi.second);
    }
    for (const auto& vi : m_actions) {
      const auto& i = vi.first;
      if (m_hermitian)
        m_action[m_index][i] = m_action[i][m_index];
      else
        m_action[m_index][i] = resres ? m_handlers->rq().dot(action, vi.second) : m_handlers->rq().dot(vector, vi.second);
    }
    m_metric[m_index][m_index] = m_handlers->rr().dot(vector, vector);
//    std::cout << "Length of new Q vector "<<m_metric[m_index][m_index]<<std::endl;
    m_action[m_index][m_index] = resres ? m_handlers->rr().dot(action, action) : m_handlers->rr().dot(vector, action);
    m_action_action[m_index][m_index] = m_handlers->rr().dot(action, action); // TODO retire this
    m_metric_pspace[m_index] = std::vector<scalar_type>(m_pspace.size());
    m_action_pspace[m_index] = std::vector<scalar_type>(m_pspace.size());
    for (auto i = 0; i < m_pspace.size(); i++) {
      m_metric_pspace[m_index][i] = m_handlers->rp().dot(vector,m_pspace[i]);
      m_action_pspace[m_index][i] = m_handlers->rp().dot(action,m_pspace[i]);
    }
    m_rhs[m_index] = std::vector<scalar_type>();
    for (const auto& rhs1 : rhs)
      m_rhs[m_index].push_back(m_handlers->rq().dot(vector, rhs1));
    m_vectors.emplace(m_index, m_handlers->qr().copy(vector));
    m_actions.emplace(m_index, m_handlers->qr().copy(action));
    m_index++;
    m_keys = keys();
  }

  /*!
   * @brief Add a new vector to the Q space formed as the normalised difference of two vectors
   * @param vector
   * @param action
   * @param oldvector
   * @param oldaction
   * @param rhs
   * @param resres If true, action matrix will be action.action instead of vector.action
   * @param orthogonalise If true, the new vector will be orthogonal to vector
   * @return The scale factor applied to make the new vector length 1
   */
  scalar_type add(const Rvector& vector, const Rvector& action, const Qvector& oldvector, const Qvector& oldaction,
                  const std::vector<Qvector>& rhs, bool resres = false, bool orthogonalise = true) {
    auto rr = m_handlers->rr().dot(vector, vector);
    typename Rvector::value_type scale_factor, diff_factor;
    if (resres) {
      rr = m_handlers->rr().dot(action, action);
      // FIXME can we avoid calculating norm of the oldaction. Was it calculated before?
      auto dd = m_handlers->qq().dot(oldaction, oldaction);
      auto rd = m_handlers->rq().dot(action, oldaction);
      scale_factor = 1 / std::sqrt(rr + dd - 2 * rd);
      diff_factor = 1;
    } else {
      // FIXME Again, can this dot be avoided
      auto dd = m_handlers->qq().dot(oldvector, oldvector);
      auto rd = m_handlers->rq().dot(vector, oldvector);
      //      std::cout << "dd-1=" << dd - 1 << ", rr-1=" << rr - 1 << ", rd-1=" << rd - 1 << std::endl;
      diff_factor = orthogonalise ? rr / rd : 1;
      auto norm = std::sqrt(std::max(rr - 2 * diff_factor * rd + diff_factor * diff_factor * dd, (decltype(rr))0));
      if (norm == 0) { // abandon ship
        return 0;
      } else {
        scale_factor = 1 / norm;
      }
    }
    //    std::cout << "Q.add difference, scale_factor=" << scale_factor << ", diff_factor=" << diff_factor
    //              << ", orthogonalise=" << orthogonalise << std::endl;
    // FIXME why are they passed as const, if they get modified?
    auto& v = const_cast<Rvector&>(vector);
    auto& a = const_cast<Rvector&>(action);
    m_handlers->rr().scal(scale_factor, v);
    m_handlers->rq().axpy(-diff_factor * scale_factor, oldvector, v);
    m_handlers->rr().scal(scale_factor, a);
    m_handlers->rq().axpy(-diff_factor * scale_factor, oldaction, a);
    auto actual_norm = std::sqrt(resres ? m_handlers->rr().dot(a, a) : m_handlers->rr().dot(v, v));
    //    std::cout << "actual_norm=" << actual_norm << std::endl;
    if (actual_norm > 1e-6 and
        std::abs(actual_norm - 1) > 1e-2) { // rescale because of numerical precision problems when vector
                                            //    \approx oldvector
      //    do not do it if the problem is severe, since then action will be inaccurate
      m_handlers->rr().scal(1. / actual_norm, v);
      m_handlers->rr().scal(1. / actual_norm, a);
      scale_factor /= actual_norm;
    }
    m_scale_factors[m_index] = scale_factor;
    m_diff_factors[m_index] = diff_factor;
    add(v, a, rhs, resres);
    //    std::cout << "new Q vector self-overlap=" << v.dot(v) << std::endl;
    m_handlers->rq().axpy(diff_factor * scale_factor, oldvector, v);
    m_handlers->rr().scal(1. / scale_factor, v);
    m_handlers->rq().axpy(diff_factor * scale_factor, oldaction, a);
    m_handlers->rr().scal(1. / scale_factor, a);
    //    std::cout << "created Q, scale_factor=" << scale_factor << ", diff_factor=" << diff_factor << std::endl;
    return scale_factor;
  }

  /*!
   * @brief Refresh stored interactions with P space. Must be called whenever the P space is changed.
   * @param workspace is used as scratch space, and its contents are undefined on exit unless the P space is null.
   */
  void refreshP(Rvector& workspace) {
    for (const auto& vi : m_vectors) {
      const auto& i = vi.first;
      m_metric_pspace[i].resize(m_pspace.size());
      m_action_pspace[i].resize(m_pspace.size());
      workspace = m_vectors[i];
      for (auto j = 0; j < m_pspace.size(); j++) {
        m_metric_pspace[i][j] = m_handlers->rp().dot(workspace, m_pspace[j]);
      }
      workspace = m_actions[i];
      for (auto j = 0; j < m_pspace.size(); j++) {
        m_action_pspace[i][j] = m_handlers->rp().dot(workspace, m_pspace[j]);
      }
    }
  }

  /*!
   * @brief Empty the P space
   */
  void clearP() {
    m_pspace.clear();
    m_metric_pspace.clear();
    m_action_pspace.clear();
  }
  /*!
   * @brief Remove a vector from the Q space
   */
  void remove(const size_t seq) {
    if (m_keys.size() <= seq)
      throw std::runtime_error("non-existent vector to erase");
    auto index = m_keys[seq];
    if (m_vectors.erase(index) != 1)
      throw std::runtime_error("non-existent vector to erase");
    if (m_actions.erase(index) != 1)
      throw std::runtime_error("non-existent vector to erase");
    m_metric.erase(index);
    m_action.erase(index);
    m_rhs.erase(index);
    for (const auto& vi : m_vectors) {
      const auto& i = vi.first;
      m_metric[i].erase(index);
      m_action[i].erase(index);
    }
    m_metric_pspace.erase(index);
    m_action_pspace.erase(index);
    m_keys = keys();
  }
};
} // namespace iterativesolver
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_Q_H_
