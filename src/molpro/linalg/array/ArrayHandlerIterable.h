#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERITERABLE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERITERABLE_H
#include "molpro/linalg/array/ArrayHandler.h"
#include <numeric>

namespace molpro {
namespace linalg {
namespace array {

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
  using ArrayHandler<AL, AR>::copySame;
  using ArrayHandler<AL, AR>::copy;
  using ArrayHandler<AL, AR>::fill;
  using ArrayHandler<AL, AR>::scal;
  using ArrayHandler<AL, AR>::axpy;
  using ArrayHandler<AL, AR>::dot;
  using typename ArrayHandler<AL, AR>::ProxyHandle;

  ArrayHandlerIterable() = default;
  ArrayHandlerIterable(const ArrayHandlerIterable<AL, AR> &) = default;

  ProxyHandle lazy_handle() override { return this->lazy_handle(*this); };

protected:
  using ArrayHandler<AL, AR>::lazy_handle;
  using ArrayHandler<AL, AR>::error;
  using ArrayHandler<AL, AR>::m_lazy_handles;
  AR copyRR(const AR &source) override { return source; };
  AL copyLL(const AL &source) override { return source; };

  template <typename T, typename S>
  T copyAny(const S &source) {
    T result{};
    std::copy(source.begin(), source.end(), std::back_inserter(result));
    return result;
  }
  AR copyLR(const AL &source) override { return copyAny<AR, AL>(source); }
  AL copyRL(const AR &source) override { return copyAny<AL, AR>(source); }

  template <typename T>
  void scalAny(value_type alpha, T &x) {
    for (auto &el : x)
      el *= alpha;
  }
  void scalL(value_type alpha, AL &x) override { scalAny(alpha, x); }
  void scalR(value_type alpha, AR &x) override { scalAny(alpha, x); }

  template <typename T>
  void fillAny(value_type alpha, T &x) {
    std::fill(x.begin(), x.end(), alpha);
  }
  void fillL(value_type alpha, AL &x) override { fillAny(alpha, x); }
  void fillR(value_type alpha, AR &x) override { fillAny(alpha, x); }

  template <typename T, typename S>
  void axpyAny(value_type alpha, const T &x, S &y) {
    if (y.size() > x.size())
      error("ArrayHandlerIterable::axpy() incompatible x and y arrays, y.size() > x.size()");
    auto ix = x.begin();
    for (auto iy = y.begin(); iy != y.end(); ++iy, ++ix)
      *iy += alpha * (*ix);
  }
  void axpyLL(value_type alpha, const AL &x, AL &y) override { axpyAny(alpha, x, y); }
  void axpyLR(value_type alpha, const AL &x, AR &y) override { axpyAny(alpha, x, y); }
  void axpyRL(value_type alpha, const AR &x, AL &y) override { axpyAny(alpha, x, y); }
  void axpyRR(value_type alpha, const AR &x, AR &y) override { axpyAny(alpha, x, y); }

  template <typename T, typename S>
  value_type dotAny(const T &x, const S &y) {
    if (x.size() > y.size())
      error("ArrayHandlerIterable::dot() incompatible x and y arrays, x.size() > y.size()");
    return std::inner_product(x.begin(), x.end(), y.begin(), (value_type)0);
  }

  value_type dotLL(const AL &x, const AL &y) { return dotAny(x, y); }
  value_type dotLR(const AL &x, const AR &y) { return dotAny(x, y); }
  value_type dotRL(const AR &x, const AL &y) { return dotAny(x, y); }
  value_type dotRR(const AR &x, const AR &y) { return dotAny(x, y); }

}; // namespace linalg

} // namespace array
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERITERABLE_H
