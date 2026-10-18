#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERDDISKDISTR_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERDDISKDISTR_H

#include <map>
#include <molpro/linalg/array/ArrayHandlerWithCopyFunc.h>
#include <molpro/linalg/array/util/gemm.h>
#include <string>

namespace molpro::linalg::array {
using molpro::linalg::array::util::gemm_inner_distr_distr;
using molpro::linalg::array::util::gemm_outer_distr_distr;

template <class AL, class AR = AL>
class ArrayHandlerDDiskDistr : public ArrayHandlerWithCopyFunc<AL, AR> {
public:
  using typename ArrayHandler<AL, AR>::value_type_L;
  using typename ArrayHandler<AL, AR>::value_type_R;
  using typename ArrayHandler<AL, AR>::value_type;
  using typename ArrayHandler<AL, AR>::value_type_abs;
  using typename ArrayHandler<AL, AR>::ProxyHandle;

  //! New arrays are created by AL::CreateTempCopy(), with any file next to that of the source, if it has one
  ArrayHandlerDDiskDistr()
      : ArrayHandlerWithCopyFunc<AL, AR>([](const AR &source) { return AL::CreateTempCopy(source); }) {}

  //! New arrays are created by AL::CreateTempCopy(), with any file in directory
  explicit ArrayHandlerDDiskDistr(const std::string &directory)
      : ArrayHandlerWithCopyFunc<AL, AR>(
            [directory](const AR &source) { return AL::CreateTempCopy(source, directory); }) {}

  //! Constructor taking a copy function, which can e.g. specify the file underlying the new array
  ArrayHandlerDDiskDistr(std::function<AL(const AR &)> copy_func)
      : ArrayHandlerWithCopyFunc<AL, AR>(std::move(copy_func)) {}

  ProxyHandle lazy_handle() override { return this->lazy_handle(*this); };

  using ArrayHandler<AL, AR>::lazy_handle;
  using ArrayHandler<AL, AR>::error;

  using ArrayHandlerWithCopyFunc<AL, AR>::copy;

  void copy(AL &x, const AR &y) override {
    this->m_counter->copy++;
    x.copy(y);
  };

  void scal(value_type alpha, AL &x) override {
    this->m_counter->scal++;
    x.scal(alpha);
  }

  void fill(value_type alpha, AL &x) override { x.fill(alpha); }

  void axpy(value_type alpha, const AR &x, AL &y) override {
    this->m_counter->axpy++;
    y.axpy(alpha, x);
  }

  value_type dot(const AL &x, const AR &y) override {
    this->m_counter->dot++;
    return x.dot(y);
  }

  void gemm_outer(const Matrix<value_type> alphas, const CVecRef<AR> &xx, const VecRef<AL> &yy) override {
    this->m_counter->gemm_outer++;
    gemm_outer_distr_distr(alphas, xx, yy);
  }

  Matrix<value_type> gemm_inner(const CVecRef<AL> &xx, const CVecRef<AR> &yy) override {
    auto prof = molpro::Profiler::single()->push("ArrayHandlerDDiskDistr::gemm_inner");
    this->m_counter->gemm_inner++;
    return gemm_inner_distr_distr(xx, yy);
  }

  std::map<size_t, value_type_abs> select_max_dot(size_t n, const AL &x, const AR &y) override {
    return x.select_max_dot(n, y);
  }

  std::map<size_t, value_type_abs> select(size_t n, const AL &x, bool max = false, bool ignore_sign = false) override {
    return x.select(n, max, ignore_sign);
  }
};

} // namespace molpro::linalg::array

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERDDISKDISTR_H
