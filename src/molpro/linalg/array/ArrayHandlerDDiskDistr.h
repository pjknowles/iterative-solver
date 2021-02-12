#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERDDISKDISTR_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERDDISKDISTR_H

#include <map>
#include <molpro/linalg/array/ArrayHandler.h>
#include <molpro/linalg/array/util/gemm.h>

using molpro::linalg::array::util::gemm_outer_distr_distr;
using molpro::linalg::array::util::gemm_inner_distr_distr;

namespace molpro::linalg::array {

template <class AL, class AR = AL>
class ArrayHandlerDDiskDistr : public ArrayHandler<AL, AR> {
public:
  using typename ArrayHandler<AL, AR>::value_type_L;
  using typename ArrayHandler<AL, AR>::value_type_R;
  using typename ArrayHandler<AL, AR>::value_type;
  using typename ArrayHandler<AL, AR>::value_type_abs;
  using typename ArrayHandler<AL, AR>::ProxyHandle;

  ArrayHandlerDDiskDistr()
      : ArrayHandler<AL, AR>(), m_copy_func([](const AR &source) { return AL::CreateTempCopy(source); }) {}

  /*!
   * @brief Constructor taking a copy function.
   *
   * Disk arrays may require extra arguments in the copy constructor to specify the underlying filename.
   * The user can pass copy_func which implements
   *
   * @param copy_func copy function
   */
  ArrayHandlerDDiskDistr(std::function<AL(const AR &)> copy_func)
      : ArrayHandler<AL, AR>(), m_copy_func(std::move(copy_func)){};

  ProxyHandle lazy_handle() override { return this->lazy_handle(*this); };

  using ArrayHandler<AL, AR>::lazy_handle;
  using ArrayHandler<AL, AR>::error;

  AL copy(const AR &source) override { return m_copy_func(source); };

  void copy(AL &x, const AR &y) override { x.copy(y); };

  void scal(value_type alpha, AL &x) override { x.scal(alpha); }

  void fill(value_type alpha, AL &x) override { x.fill(alpha); }

  void axpy(value_type alpha, const AR &x, AL &y) override { y.axpy(alpha, x); }

  value_type dot(const AL &x, const AR &y) override { return x.dot(y); }
  
  void gemm_outer(const Matrix<value_type> alphas, const CVecRef<AR> &xx, const VecRef<AL> &yy) override {
    gemm_outer_distr_distr(alphas, xx, yy);
  }
  
  Matrix<value_type> gemm_inner(const CVecRef<AL> &xx, const CVecRef<AR> &yy) override {
    return gemm_inner_distr_distr(xx, yy);
  }
  
  std::map<size_t, value_type_abs> select_max_dot(size_t n, const AL &x, const AR &y) override {
    return x.select_max_dot(n, y);
  }

protected:
  std::function<AL(const AR &)> m_copy_func; //!< function for making a copy
};

} // namespace molpro::linalg::array

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERDDISKDISTR_H
