#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_P_H_
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_P_H_
#include <map>
#include <molpro/linalg/array/ArrayHandler.h>

namespace molpro {
namespace linalg {
namespace itsolv {

template <class Pvector = std::map<size_t, double>>
class P {
public:
  using value_type = typename Pvector::mapped_type;

protected:
  std::vector<value_type> m_metric;
  std::vector<value_type> m_action;
  std::vector<value_type> m_rhs;
  std::vector<Pvector> m_vectors;

public:
  P() {}

  const value_type& metric(size_t i, size_t j) const { return m_metric[m_vectors.size() * j + i]; }

  const value_type& action(size_t i, size_t j) const { return m_action[m_vectors.size() * j + i]; }

  const value_type& rhs(size_t i, size_t j) const { return m_rhs[m_vectors.size() * j + i]; }

  size_t size() const { return m_vectors.size(); }

  const Pvector& operator[](size_t i) const { return m_vectors[i]; }

  /*!
   * @brief Add a new vector to the space. Also compute and store the new elements of the overlap and action
   * matrices
   * @param Pvectors the vectors to add
   * @param PP Matrix projected onto the existing+new, new P space. It should be provided as a
   * 1-dimensional array, with the existing+new index running fastest.
   * @param rhs The right hand side of an inhomogeneous equation set. The projection onto the new P vectors will be
   * calculated.
   * @param handle_pp handle between arrays in P space
   * @param handle_qp handle between arrays Q and P space
   * @tparam Qvector Out-of memory vectors to be used for storing the space.
   * @tparam Pvector Containers that specify a P-space member
   */
  template <class Qvector>
  void add(const std::vector<Pvector>& Pvectors, const value_type* PP, const std::vector<Qvector>& rhs,
           array::ArrayHandler<Pvector, Pvector>& handle_pp, array::ArrayHandler<Qvector, Pvector>& handle_qp) {
    auto old_size = m_vectors.size();
    auto new_size = m_vectors.size() + Pvectors.size();
    {
      std::vector<value_type> new_metric(new_size * new_size);
      std::vector<value_type> new_action(new_size * new_size);
      std::vector<value_type> new_rhs(new_size * rhs.size());
      for (size_t i = 0; i < old_size; i++) {
        for (size_t j = 0; j < old_size; j++) {
          new_metric[i * new_size + j] = m_metric[i * old_size + j];
          new_action[i * new_size + j] = m_action[i * old_size + j];
        }
        for (size_t j = 0; j < rhs.size(); j++) {
          new_rhs[i + new_size * j] = m_action[i + old_size * j];
        }
      }
      for (size_t i = 0; i < Pvectors.size(); i++) {
        for (size_t j = 0; j < new_size; j++)
          new_action[j * new_size + (i + old_size)] = new_action[j + new_size * (i + old_size)] = PP[new_size * j + i];
        m_vectors.emplace_back(handle_pp.copy(Pvectors[i]));
        for (size_t j = 0; j < m_vectors.size(); j++) {
          auto overlap = handle_pp.dot(Pvectors[i], m_vectors[j]);
          new_metric[j * new_size + (i + old_size)] = new_metric[j + new_size * (i + old_size)] = overlap;
        }
        for (size_t j = 0; j < rhs.size(); j++) {
          new_rhs[i + old_size + new_size * j] = handle_qp.dot(rhs[j], m_vectors[i]);
        }
      }
      m_metric = new_metric;
      m_action = new_action;
      m_rhs = new_rhs;
    }
  }
  void clear() {
    m_metric.clear();
    m_action.clear();
    m_rhs.clear();
    m_vectors.clear();
  }
};
} // namespace itsolv
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_P_H_
