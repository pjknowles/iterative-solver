#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERITERABLE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERITERABLE_H
#include <cstddef>
#include <molpro/linalg/array/ArrayHandler.h>
#include <molpro/linalg/array/util.h>
#include <molpro/linalg/array/util/gemm.h>
#include <molpro/linalg/array/util/select.h>
#include <molpro/linalg/array/util/select_max_dot.h>
#include <numeric>
#include <stdexcept>
#include <string>

namespace molpro::linalg::array {
using molpro::linalg::array::util::gemm_inner_default;
using molpro::linalg::array::util::gemm_outer_default;

/*!
 * @brief Array handler for two containers both of which can be iterated through using begin() and end() member
 * functions, and have a copy constructor.
 */
template <typename AL, typename AR = AL>
class ArrayHandlerIterable : public ArrayHandler<AL, AR> {
public:
  using typename ArrayHandler<AL, AR>::value_type_L;
  using typename ArrayHandler<AL, AR>::value_type_R;
  using typename ArrayHandler<AL, AR>::value_type;
  using typename ArrayHandler<AL, AR>::value_type_abs;
  using typename ArrayHandler<AL, AR>::ProxyHandle;

  ArrayHandlerIterable() = default;
  ArrayHandlerIterable(const ArrayHandlerIterable<AL, AR> &) = default;

  AL copy(const AR &source) override {
    auto result = util::allocate_like<AL>(source.size());
    copy(result, source);
    return result;
  };

  void copy(AL &x, const AR &y) override {
    using std::begin;
    using std::end;
    std::copy(begin(y), end(y), begin(x));
  };

  void scal(value_type alpha, AL &x) override {
    for (auto &el : x)
      el *= alpha;
  };

  void fill(value_type alpha, AL &x) override {
    using std::begin;
    using std::end;
    std::fill(begin(x), end(x), alpha);
  };

  void axpy(value_type alpha, const AR &x, AL &y) override {
    auto prof = molpro::Profiler::single();
    prof->start("ArrayHandlerIterable::axpy");
    if (x.size() < y.size())
      error("ArrayHandlerIterable::axpy() incompatible x and y arrays, x.size() < y.size()");
    using std::begin;
    using std::end;
    std::transform(begin(y), end(y), begin(x), begin(y), [alpha](auto &ely, auto &elx) { return ely + alpha * elx; });
    prof->stop();
  };

  /*!
   * @brief The hermitian inner product <x|y>, i.e. conjugate-linear in x and linear in y.
   *
   * The real case is selected at compile time and is the plain product this function computed before
   * conjugation was introduced. A release build optimises the conjugating lambda away entirely -- the
   * two are within noise of each other at -O3 on both GCC and clang -- but at -O0 nothing is inlined
   * and the extra call per element costs slightly over a factor of two on this loop, so the branch is
   * worth the four lines.
   */
  value_type dot(const AL &x, const AR &y) override {
    if (x.size() > y.size())
      error("ArrayHandlerIterable::dot() incompatible x and y arrays, x.size() > y.size()");
    using std::begin;
    using std::end;
    if constexpr (molpro::linalg::is_complex<value_type>{})
      return std::inner_product(begin(x), end(x), begin(y), value_type{}, std::plus<value_type>{},
                                [](const auto &elx, const auto &ely) {
                                  return molpro::linalg::conjugate(static_cast<value_type>(elx)) *
                                         static_cast<value_type>(ely);
                                });
    else
      return std::inner_product(begin(x), end(x), begin(y), value_type{});
  };

  void gemm_outer(const Matrix<value_type> alphas, const CVecRef<AR> &xx, const VecRef<AL> &yy) override {
    gemm_outer_default(*this, alphas, xx, yy);
  }

  Matrix<value_type> gemm_inner(const CVecRef<AL> &xx, const CVecRef<AR> &yy) override {
    return gemm_inner_default(*this, xx, yy);
  }

  std::map<size_t, value_type_abs> select_max_dot(size_t n, const AL &x, const AR &y) override {
    if (n > x.size() || n > y.size())
      error("ArrayHandlerIterable::select_max_dot() n is too large");
    return util::select_max_dot<AL, AR, value_type, value_type_abs>(n, x, y);
  }

  std::map<size_t, value_type> select(size_t n, const AL &x, bool max = false, bool ignore_sign = false) override {
    if (n > x.size())
      error("ArrayHandlerIterable::select() n is too large");
    return util::select<AL, value_type>(n, x, max, ignore_sign);
  }

  ProxyHandle lazy_handle() override { return this->lazy_handle(*this); };

protected:
  using ArrayHandler<AL, AR>::error;
  using ArrayHandler<AL, AR>::lazy_handle;
  using ArrayHandler<AL, AR>::m_lazy_handles;
};

} // namespace molpro::linalg::array

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERITERABLE_H
