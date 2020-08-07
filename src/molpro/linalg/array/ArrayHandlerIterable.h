#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERITERABLE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERITERABLE_H
#include "molpro/linalg/array/ArrayHandler.h"
#include <numeric>

namespace molpro {
namespace linalg {
namespace array {
namespace util {

template <typename T, bool = std::is_array<T>::value>
struct is_array : std::true_type {};

template <typename T>
struct is_array<T, false> : std::false_type {};

template <typename T, int N>
struct is_array<std::array<T, N>, false> : std::true_type {};

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
  using typename ArrayHandler<AL, AR>::ProxyHandle;

  ArrayHandlerIterable() = default;
  ArrayHandlerIterable(const ArrayHandlerIterable<AL, AR> &) = default;

  AL copy(const AR &source) override { return copyAny<AL, AR>(source); };

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

  ProxyHandle lazy_handle() override { return this->lazy_handle(*this); };

protected:
  using ArrayHandler<AL, AR>::error;
  using ArrayHandler<AL, AR>::lazy_handle;
  using ArrayHandler<AL, AR>::m_lazy_handles;

  template <typename T, typename S, typename std::enable_if_t<util::is_array<T>::value, nullptr_t> = nullptr>
  T copyAny(const S &source) {
    auto result = T();
    std::copy(begin(source), end(source), begin(result));
    return result;
  }
  template <typename T, typename S, typename std::enable_if_t<!util::is_array<T>::value, int> = 0>
  T copyAny(const S &source) {
    auto result = T(source.size());
    std::copy(begin(source), end(source), begin(result));
    return result;
  }
}; // namespace linalg

} // namespace array
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERITERABLE_H
