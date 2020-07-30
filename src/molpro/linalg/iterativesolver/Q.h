#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_Q_H_
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_Q_H_
#include "P.h"
#include <cmath>
#include <map>
#include <memory>
#include <molpro/iostream.h>
#include <vector>

/*!
 * @brief A class holding a Q space
 * @tparam fastvector In-memory vectors. New instances will not be created.
 * @tparam slowvector Out-of memory vectors to be used for storing the space.
 */
template <class fastvector, class slowvector = fastvector>
class Q {
  using scalar_type = decltype(
      std::declval<fastvector>().dot(std::declval<const fastvector&>())); ///< The type of scalar products of vectors
  using key_type = int;
  bool m_hermitian;
  std::map<key_type, std::map<key_type, scalar_type>> m_metric;
  std::map<key_type, std::map<key_type, scalar_type>> m_action;
  std::map<key_type, std::map<key_type, scalar_type>> m_action_action;
  std::map<key_type, std::vector<scalar_type>> m_rhs;
  key_type m_index = 0;
  std::map<key_type, slowvector> m_vectors;
  std::map<key_type, slowvector> m_actions;
  const P<typename fastvector::value_type>& m_pspace;
  std::map<key_type, std::vector<scalar_type>> m_metric_pspace;
  std::map<key_type, std::vector<scalar_type>> m_action_pspace;
  std::vector<key_type> m_keys;
  std::map<key_type, typename fastvector::value_type> m_scale_factors, m_diff_factors;

public:
  Q(const P<typename fastvector::value_type>& pspace, bool hermitian = false)
      : m_hermitian(hermitian), m_pspace(pspace) {}

  const scalar_type& metric(size_t i, size_t j) const { return m_metric.at(m_keys[i]).at(m_keys[j]); }

  const scalar_type& action(size_t i, size_t j) const { return m_action.at(m_keys[i]).at(m_keys[j]); }
  const scalar_type& action_action(size_t i, size_t j) const { return m_action_action.at(m_keys[i]).at(m_keys[j]); }

  const std::vector<scalar_type>& rhs(size_t i) const { return m_rhs.at(m_keys[i]); }

  const std::vector<scalar_type>& metric_pspace(size_t i) const { return m_metric_pspace.at(m_keys[i]); }
  const std::vector<scalar_type>& action_pspace(size_t i) const { return m_action_pspace.at(m_keys[i]); }

  size_t size() const { return m_vectors.size(); }

  const slowvector& operator[](size_t i) const { return m_vectors.at(m_keys[i]); }
  const slowvector& action(size_t i) const { return m_actions.at(m_keys[i]); }

  const typename fastvector::value_type scale_factor(size_t i) const { return m_scale_factors.at(m_keys[i]); }

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
  void add(const fastvector& vector, const fastvector& action, const std::vector<slowvector>& rhs,
           bool resres = false) {
    for (const auto& vi : (resres ? m_actions : m_vectors)) {
      m_metric[m_index][vi.first] = m_metric[vi.first][m_index] = vector.dot(vi.second);
      m_action[vi.first][m_index] = action.dot(vi.second);
    }
    for (const auto& vi : m_actions) {
      m_action_action[m_index][vi.first] = m_action_action[vi.first][m_index] = action.dot(vi.second);
    }
    for (const auto& vi : m_actions) {
      const auto& i = vi.first;
      if (m_hermitian)
        m_action[m_index][i] = m_action[i][m_index];
      else
        m_action[m_index][i] = resres ? action.dot(vi.second) : vector.dot(vi.second);
    }
    m_metric[m_index][m_index] = vector.dot(vector);
    m_action[m_index][m_index] = resres ? action.dot(action) : vector.dot(action);
    m_action_action[m_index][m_index] = action.dot(action); // TODO retire this
    m_metric_pspace[m_index] = std::vector<scalar_type>(m_pspace.size());
    m_action_pspace[m_index] = std::vector<scalar_type>(m_pspace.size());
    for (auto i = 0; i < m_pspace.size(); i++) {
      m_metric_pspace[m_index][i] = vector.dot(m_pspace[i]);
      m_action_pspace[m_index][i] = action.dot(m_pspace[i]);
    }
    m_rhs[m_index] = std::vector<scalar_type>();
    for (const auto& rhs1 : rhs)
      m_rhs[m_index].push_back(vector.dot(rhs1));
    m_vectors.emplace(std::make_pair(m_index, slowvector{vector}));
    m_actions.emplace(std::make_pair(m_index, slowvector{action}));
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
  scalar_type add(const fastvector& vector, const fastvector& action, const slowvector& oldvector,
                  const slowvector& oldaction, const std::vector<slowvector>& rhs, bool resres = false,
                  bool orthogonalise = true) {
    auto rr = vector.dot(vector);
    typename fastvector::value_type scale_factor, diff_factor;
    if (resres) {
      rr = action.dot(action);
      auto dd = oldaction.dot(oldaction);
      auto rd = action.dot(oldaction);
      scale_factor = 1 / std::sqrt(rr + dd - 2 * rd);
      diff_factor = scale_factor;
    } else {
      auto dd = oldvector.dot(oldvector);
      auto rd = vector.dot(oldvector);
      //      molpro::cout << "dd="<<dd<<std::endl;
      //      molpro::cout << "rd="<<rd<<std::endl;
      //      molpro::cout << "rr="<<rr<<std::endl;
      diff_factor = orthogonalise ? rr / rd : 1;
      auto norm = std::sqrt(std::max(rr - 2 * diff_factor * rd + diff_factor * diff_factor * dd, (decltype(rr))0));
      if (norm == 0) { // let linear dependence code deal with this exceptional case later
        scale_factor = 1;
      } else {
        scale_factor = 1 / norm;
      }
    }
    //    molpro::cout << "Q.add difference, alpha=" << alpha << ", beta=" << beta << std::endl;
    //    molpro::cout << "dd=" << dd << ", rr=" << rr
    //                 << ", rd=" << rd << std::endl;
    auto& v = const_cast<fastvector&>(vector);
    auto& a = const_cast<fastvector&>(action);
    v.scal(scale_factor);
    v.axpy(-diff_factor * scale_factor, oldvector);
    a.scal(scale_factor);
    a.axpy(-diff_factor * scale_factor, oldaction);
    m_scale_factors[m_index] = scale_factor;
    m_diff_factors[m_index] = diff_factor;
    add(v, a, rhs, resres);
    v.axpy(diff_factor * scale_factor, oldvector);
    v.scal(1 / scale_factor);
    a.axpy(diff_factor * scale_factor, oldaction);
    a.scal(1 / scale_factor);
    //    molpro::cout << "created Q, scale_factor="<<scale_factor<<std::endl;
    return scale_factor;
  }

  /*!
   * @brief Refresh stored interactions with P space. Must be called whenever the P space is changed.
   * @param workspace is used as scratch space, and its contents are undefined on exit unless the P space is null.
   */
  void refreshP(fastvector& workspace) {
    for (const auto& vi : m_vectors) {
      const auto& i = vi.first;
      m_metric_pspace[i].resize(m_pspace.size());
      m_action_pspace[i].resize(m_pspace.size());
      workspace = m_vectors[i];
      for (auto j = 0; j < m_pspace.size(); j++) {
        m_metric_pspace[i][j] = workspace.dot(m_pspace[j]);
      }
      workspace = m_actions[i];
      for (auto j = 0; j < m_pspace.size(); j++) {
        m_action_pspace[i][j] = workspace.dot(m_pspace[j]);
      }
    }
  }

  /*!
   * @brief Remove a vector from the Q space
   * @param index
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

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_Q_H_
