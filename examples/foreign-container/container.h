#ifndef LINEARALGEBRA_CONTAINER_H
#define LINEARALGEBRA_CONTAINER_H
#include <algorithm>
#include <map>
#include <molpro/linalg/array/util/select.h>
#include <numeric>
#include <stdexcept>
#include <vector>
#include <iostream>

template <typename T = double>
class container {
  std::vector<T> m_data;

public:
  container(size_t n) : m_data(n) {}

  container(const std::map<size_t, T> &source) { throw std::logic_error("unimplementable"); }

  container& operator=(const std::map<size_t, T> &source) {
    m_data.assign(m_data.size(), 0);
    for (const auto &s : source)
      m_data[s.first] = s.second;
    return *this;
  }

  using value_type = typename decltype(m_data)::value_type;

  std::vector<value_type> &data() { return m_data; }

  const std::vector<value_type> &data() const { return m_data; }

  void fill(value_type a) {
    std::for_each(m_data.begin(), m_data.end(), [a](value_type &y) { y = a; });
  }

  void scal(value_type a) {
    std::for_each(m_data.begin(), m_data.end(), [a](value_type &y) { y *= a; });
  }

  void axpy(value_type a, const container &x) {
    std::transform(x.m_data.begin(), x.m_data.end(), m_data.begin(), m_data.begin(),
                   [a](const value_type &xx, const value_type &yy) { return yy + a * xx; });
  }

  value_type dot(const container &x) const {
    return std::inner_product(m_data.begin(), m_data.end(), x.m_data.begin(), static_cast<value_type>(0));
  }

  std::map<size_t, value_type> select_max_dot(size_t n, const container &y) const {
    throw std::logic_error("unimplemented");
  }

  std::map<size_t, value_type> select(size_t n, bool max = false, bool ignore_sign = false) const {
    return molpro::linalg::array::util::select<decltype(m_data), value_type>(n, m_data, max, ignore_sign);
  }
};

#endif // LINEARALGEBRA_CONTAINER_H
