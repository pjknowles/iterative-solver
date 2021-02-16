#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERDISTR_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERDISTR_H

#include <molpro/linalg/array/ArrayHandler.h>
#include <molpro/linalg/array/util/gemm.h>

using molpro::linalg::array::util::gemm_outer_distr_distr;
using molpro::linalg::array::util::gemm_inner_distr_distr;

namespace molpro::linalg::array {

template <class AL, class AR = AL>
class ArrayHandlerDistr : public ArrayHandler<AL, AR> {
public:
  using typename ArrayHandler<AL, AR>::value_type_L;
  using typename ArrayHandler<AL, AR>::value_type_R;
  using typename ArrayHandler<AL, AR>::value_type;
  using typename ArrayHandler<AL, AR>::value_type_abs;
  using typename ArrayHandler<AL, AR>::ProxyHandle;

  ProxyHandle lazy_handle() override { return this->lazy_handle(*this); };

  using ArrayHandler<AL, AR>::lazy_handle;
  using ArrayHandler<AL, AR>::error;
  using ArrayHandler<AL, AR>::counter;

  AL copy(const AR &source) override { this->m_counter->copy++; return AL{source}; };

  void copy(AL &x, const AR &y) override { this->m_counter->copy++; x.copy(y); };

  void scal(value_type alpha, AL &x) override { this->m_counter->scal++; x.scal(alpha); }

  void fill(value_type alpha, AL &x) override { x.fill(alpha); }

  void axpy(value_type alpha, const AR &x, AL &y) override { this->m_counter->axpy++; y.axpy(alpha, x); }

  value_type dot(const AL &x, const AR &y) override { this->m_counter->dot++; return x.dot(y); }
  
  void gemm_outer(const Matrix<value_type> alphas, const CVecRef<AR> &xx, const VecRef<AL> &yy) override {
    this->m_counter->gemm_outer++;
    gemm_outer_distr_distr(alphas, xx, yy);
  }
  
  Matrix<value_type> gemm_inner(const CVecRef<AL> &xx, const CVecRef<AR> &yy) override {
    this->m_counter->gemm_inner++;
    return gemm_inner_distr_distr(xx, yy);
  }
  
  std::map<size_t, value_type_abs> select_max_dot(size_t n, const AL &x, const AR &y) override {
    return x.select_max_dot(n, y);
  }
};

} // namespace molpro::linalg::array
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERDISTR_H
