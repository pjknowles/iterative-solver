#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERITERABLE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERITERABLE_H
#include <molpro/linalg/array/ArrayHandler.h>
#include <molpro/linalg/array/util/select_max_dot.h>
#include <molpro/linalg/array/util/gemm.h>
#include <cstddef>
#include <numeric>

using molpro::linalg::array::util::gemm_outer_default;
using molpro::linalg::array::util::gemm_inner_default;

namespace molpro::linalg::array {
namespace util {

template <typename T>
struct is_std_array : std::false_type {};

template <typename T, std::size_t N>
struct is_std_array<std::array<T, N>> : std::true_type {};

template <typename T>
constexpr bool is_std_array_v = is_std_array<T>::value;

template <typename T>
constexpr bool is_array_v = std::is_array<T>::value || is_std_array_v<T>;

} // namespace util

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

  AL copy(const AR &source) override { return copyAny<AL, AR>(source); };

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
    if (x.size() < y.size())
      error("ArrayHandlerIterable::axpy() incompatible x and y arrays, x.size() < y.size()");
    using std::begin;
    using std::end;
    std::transform(begin(y), end(y), begin(x), begin(y), [alpha](auto &ely, auto &elx) { return ely + alpha * elx; });
  };

  value_type dot(const AL &x, const AR &y) override {
    if (x.size() > y.size())
      error("ArrayHandlerIterable::dot() incompatible x and y arrays, x.size() > y.size()");
    using std::begin;
    using std::end;
    return std::inner_product(begin(x), end(x), begin(y), (value_type)0);
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

  ProxyHandle lazy_handle() override { return this->lazy_handle(*this); };

protected:
  using ArrayHandler<AL, AR>::error;
  using ArrayHandler<AL, AR>::lazy_handle;
  using ArrayHandler<AL, AR>::m_lazy_handles;

  template <typename T, typename S, typename std::enable_if_t<util::is_array_v<T>, std::nullptr_t> = nullptr>
  T copyAny(const S &source) {
    auto result = T();
    using std::begin;
    using std::end;
    std::copy(begin(source), end(source), begin(result));
    return result;
  }

  template <typename T, typename S, typename std::enable_if_t<!util::is_array_v<T>, int> = 0>
  T copyAny(const S &source) {
    auto result = T(source.size());
    using std::begin;
    using std::end;
    std::copy(begin(source), end(source), begin(result));
    return result;
  }
};

} // namespace molpro::linalg::array

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERITERABLE_H
