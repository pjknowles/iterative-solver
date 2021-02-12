#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERITERABLESPARSE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERITERABLESPARSE_H
#include <molpro/linalg/array/ArrayHandler.h>
#include <molpro/linalg/array/util/select_max_dot.h>
#include <molpro/linalg/array/util/gemm.h>
#include <numeric>

using molpro::linalg::array::util::gemm_outer_default;
using molpro::linalg::array::util::gemm_inner_default;

namespace molpro::linalg::array {
template <typename AL, typename AR, bool = has_mapped_type_v<AR>>
class ArrayHandlerIterableSparse : public ArrayHandler<AL, AR> {};

/*!
 * @brief Array handler between an iterable and a sparse arrays. Iterable container must implement access operator.
 */
template <typename AL, typename AR>
class ArrayHandlerIterableSparse<AL, AR, true> : public ArrayHandler<AL, AR> {
public:
  using typename ArrayHandler<AL, AR>::value_type_L;
  using typename ArrayHandler<AL, AR>::value_type_R;
  using typename ArrayHandler<AL, AR>::value_type;
  using typename ArrayHandler<AL, AR>::value_type_abs;
  using typename ArrayHandler<AL, AR>::ProxyHandle;

  AL copy(const AR &source) override {
    static_assert(true, "General copy from sparse to dense is ill-defined");
    return AL{};
  };

  void copy(AL &x, const AR &y) override { static_assert(true, "General copy from sparse to dense is ill-defined"); }

  void scal(value_type alpha, AL &x) override { static_assert(true, "Use ArrayHandlerIterable for unary operations"); };

  void fill(value_type alpha, AL &x) override { static_assert(true, "Use ArrayHandlerIterable for unary operations"); };

  void axpy(value_type alpha, const AR &x, AL &y) override {
    for (const auto &el_x : x)
      if (el_x.first < y.size())
        y[el_x.first] += alpha * el_x.second;
  };

  value_type dot(const AL &x, const AR &y) override {
    value_type tot = 0;
    for (const auto &el_y : y)
      if (el_y.first < x.size())
        tot += x[el_y.first] * el_y.second;
    return tot;
  };
  
  void gemm_outer(const Matrix<value_type> alphas, const CVecRef<AR> &xx, const VecRef<AL> &yy) override {
    gemm_outer_default(*this, alphas, xx, yy);
  }
  
  Matrix<value_type> gemm_inner(const CVecRef<AL> &xx, const CVecRef<AR> &yy) override {
    return gemm_inner_default(*this, xx, yy);
  }
  
  std::map<size_t, value_type_abs> select_max_dot(size_t n, const AL &x, const AR &y) override {
    if (n > x.size() || n > y.size())
      error("ArrayHandlerIterableSparse::select_max_dot() n is too large");
    return util::select_max_dot_iter_sparse<AL, AR, value_type, value_type_abs>(n, x, y);
  }

  ProxyHandle lazy_handle() override { return this->lazy_handle(*this); };

protected:
  using ArrayHandler<AL, AR>::error;
  using ArrayHandler<AL, AR>::lazy_handle;
  using ArrayHandler<AL, AR>::m_lazy_handles;
};

} // namespace molpro::linalg::array

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERITERABLESPARSE_H
