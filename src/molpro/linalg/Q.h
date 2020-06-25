#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_Q_H_
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_Q_H_
#include "P.h"
#include <cmath>
#include <map>
#include <memory>

template <class fastvector, class slowvector = fastvector>
class Q {
  using scalar_type = decltype(
      std::declval<fastvector>().dot(std::declval<const fastvector&>())); ///< The type of scalar products of vectors
  bool m_hermitian;
  std::map<int, std::map<int, scalar_type>> m_metric;
  std::map<int, std::map<int, scalar_type>> m_action;
  int m_index = 0;
  std::map<int, slowvector> m_vectors;
  std::map<int, slowvector> m_actions;

public:
  std::shared_ptr<P<typename fastvector::value_type, scalar_type>> m_pspace;

public:
  Q(std::shared_ptr<P<typename fastvector::value_type, scalar_type>> pspace = nullptr, bool hermitian = false)
      : m_hermitian(hermitian or pspace != nullptr), m_pspace(pspace) {}

  scalar_type metric(int i, int j) const { return m_metric[i][j]; }

  scalar_type action(int i, int j) const { return m_action[i][j]; }

  void add(const std::vector<fastvector>& vectors, const std::vector<fastvector>& actions) {
    assert(vectors.size() == actions.size());
    for (auto i = 0; i < vectors.size(); i++)
      add(vectors[i], actions[i]);
  }

  /*!
   * @brief Add a new vector to the Q space. Also compute and store the new elements of the QQ overlap and action
   * matrices
   * @param vector
   * @param action
   */
  void add(const fastvector& vector, const fastvector& action) {
    for (const auto& vi : m_vectors) {
      const auto& i = vi.first;
      m_metric[m_index][i] = m_metric[i][m_index] = vector.dot(vi.second);
      m_action[i][m_index] = action.dot(vi.second);
    }
    for (const auto& vi : m_actions) {
      const auto& i = vi.first;
      if (m_hermitian)
        m_action[m_index][i] = m_action[i][m_index];
      else
        m_action[m_index][i] = vector.dot(vi.second);
    }
    m_metric[m_index][m_index] = vector.dot(vector);
    m_action[m_index][m_index] = vector.dot(action);
    m_vectors.emplace(std::make_pair(m_index, slowvector{vector}));
    m_actions.emplace(std::make_pair(m_index, slowvector{action}));
    m_index++;
  }

  /*!
   * @brief Add a new vector to the Q space formed as the normalised difference of two vectors
   * @param vector
   * @param action
   * @param oldvector
   * @param oldaction
   */
  void add(const fastvector& vector, const fastvector& action, const slowvector& oldvector,
           const slowvector& oldaction) {
    auto norm = vector.dot(vector) + oldvector.dot(oldvector) - 2 * vector.dot(oldvector);
    auto scale = 1 / std::sqrt(norm);
    auto& v = const_cast<fastvector>(vector);
    auto& a = const_cast<fastvector>(action);
    v.scal(scale);
    v.axpy(-scale, oldvector);
    a.scal(scale);
    a.axpy(-scale, oldaction);
    add(v, a);
    v.axpy(scale, oldvector);
    v.scal(1 / scale);
    a.axpy(scale, oldaction);
    a.scal(1 / scale);
  }

  /*!
   * @brief Remove a vector from the Q space
   * @param index
   */
  void remove(const int index) {
    if (m_vectors.erase(index) != 1)
      throw std::runtime_error("non-existent vector to erase");
    if (m_actions.erase(index) != 1)
      throw std::runtime_error("non-existent vector to erase");
    m_metric.erase(index);
    m_action.erase(index);
    for (const auto& vi : m_vectors) {
      const auto& i = vi.first;
      m_metric[i].erase(index);
      m_action[i].erase(index);
    }
  }
};

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_Q_H_
