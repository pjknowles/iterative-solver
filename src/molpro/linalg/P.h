#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_P_H_
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_P_H_
#include <map>

template <class value_type = double, class scalar_type = value_type>
class P {
  using Pvector = std::map<size_t, scalar_type>;
//  std::map<int, std::map<int, scalar_type>> m_metric;
//  std::map<int, std::map<int, scalar_type>> m_action;
  std::vector<scalar_type> m_metric;
  std::vector<scalar_type> m_action;
  std::map<int, Pvector> m_vectors;

public:
  P() :  {}

  scalar_type metric(int i, int j) const { return m_metric[m_vectors.size()*j+i]; }

  scalar_type action(int i, int j) const { return m_action[m_vectors.size()*j+i]; }

  void add(const std::vector<Pvector>& vectors, const std::vector<Pvector>& actions) {
    assert(vectors.size() == actions.size());
    for (auto i = 0; i < vectors.size(); i++)
      add(vectors[i], actions[i]);
  }

  /*!
   * @brief Add a new vector to the space. Also compute and store the new elements of the overlap and action
   * matrices
   * @param Pvectors the vectors to add
   * @param PP Matrix projected onto the existing+new, new P space. It should be provided as a
   */
  void add(const std::vector<Pvector>& Pvectors, const scalar_type* PP) {
    auto new_size = m_vectors.size() + Pvectors.size();
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
   * @brief Remove a vector from the space
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

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_P_H_
