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
template <typename AL, typename AR = AL> class ArrayHandlerIterable : public ArrayHandler<AL, AR> {
public:
  using typename ArrayHandler<AL, AR>::value_type;
  using ArrayHandler<AL, AR>::copySame;
  using ArrayHandler<AL, AR>::copy;
  using ArrayHandler<AL, AR>::fill;
  using ArrayHandler<AL, AR>::scal;
  using ArrayHandler<AL, AR>::axpy;
  using ArrayHandler<AL, AR>::dot;

  ArrayHandlerIterable() = default;
  ArrayHandlerIterable(const ArrayHandlerIterable<AL, AR> &) = default;

  void error(std::string message) override { throw std::runtime_error(message); }

  typename ArrayHandler<AL, AR>::ProxyHandle lazy_handle() override {
    auto handle = std::make_shared<typename ArrayHandler<AL, AR>::LazyHandle>(*this);
    save_handle(handle);
    return handle;
  };

protected:
  using ArrayHandler<AL, AR>::save_handle;
  using ArrayHandler<AL, AR>::m_lazy_handles;
  AR copyRR(const AR &source) override { return source; };
  AL copyLL(const AL &source) override { return source; };

  template <typename T, typename S> T copyAny(const S &source) {
    T result{};
    std::copy(source.begin(), source.end(), std::back_inserter(result));
    return result;
  }
  AR copyLR(const AL &source) override { return copyAny<AR, AL>(source); }
  AL copyRL(const AR &source) override { return copyAny<AL, AR>(source); }

  template <typename T> T &scalAny(value_type alpha, T &x) {
    for (auto &el : x)
      el *= alpha;
    return x;
  }
  AL &scalL(value_type alpha, AL &x) override { return scalAny(alpha, x); }
  AR &scalR(value_type alpha, AR &x) override { return scalAny(alpha, x); }

  template <typename T> T &fillAny(value_type alpha, T &x) {
    std::fill(x.begin(), x.end(), alpha);
    return x;
  }
  AL &fillL(value_type alpha, AL &x) override { return fillAny(alpha, x); }
  AR &fillR(value_type alpha, AR &x) override { return fillAny(alpha, x); }

  template <typename T, typename S> S &axpyAny(value_type alpha, const T &x, S &y) {
    if (y.size() > x.size())
      error("ArrayHandlerIterable::axpy() incompatible x and y arrays, y.size() > x.size()");
    auto ix = x.begin();
    for (auto iy = y.begin(); iy != y.end(); ++iy, ++ix)
      *iy += alpha * (*ix);
    return y;
  }
  AL &axpyLL(value_type alpha, const AL &x, AL &y) override { return axpyAny(alpha, x, y); }
  AR &axpyLR(value_type alpha, const AL &x, AR &y) override { return axpyAny(alpha, x, y); }
  AL &axpyRL(value_type alpha, const AR &x, AL &y) override { return axpyAny(alpha, x, y); }
  AR &axpyRR(value_type alpha, const AR &x, AR &y) override { return axpyAny(alpha, x, y); }

  template <typename T, typename S> value_type dotAny(const T &x, const S &y) {
    if (x.size() > y.size())
      error("ArrayHandlerIterable::dot() incompatible x and y arrays, x.size() > y.size()");
    return std::inner_product(x.begin(), x.end(), y.begin(), (value_type)0);
  }

  value_type dotLL(const AL &x, const AL &y) { return dotAny(x, y); }
  value_type dotLR(const AL &x, const AR &y) { return dotAny(x, y); }
  value_type dotRL(const AR &x, const AL &y) { return dotAny(x, y); }
  value_type dotRR(const AR &x, const AR &y) { return dotAny(x, y); }

  template <typename T, typename S>
  void fused_axpyAny(const std::vector<std::tuple<size_t, size_t, size_t>> &reg, const std::vector<value_type> &alphas,
                     const std::vector<std::reference_wrapper<const T>> &xx,
                     std::vector<std::reference_wrapper<S>> &yy) {
    for (size_t i = 0; i < reg.size(); ++i) {
      size_t ai, xi, yi;
      std::tie(ai, xi, yi) = reg[i];
      axpy(alphas[ai], xx[xi].get(), yy[yi].get());
    }
  }
  void fused_axpyLL(const std::vector<std::tuple<size_t, size_t, size_t>> &reg, const std::vector<value_type> &alphas,
                    const std::vector<std::reference_wrapper<const AL>> &xx,
                    std::vector<std::reference_wrapper<AL>> &yy) override {
    return fused_axpyAny(reg, alphas, xx, yy);
  }
  void fused_axpyLR(const std::vector<std::tuple<size_t, size_t, size_t>> &reg, const std::vector<value_type> &alphas,
                    const std::vector<std::reference_wrapper<const AL>> &xx,
                    std::vector<std::reference_wrapper<AR>> &yy) override {
    return fused_axpyAny(reg, alphas, xx, yy);
  }
  void fused_axpyRL(const std::vector<std::tuple<size_t, size_t, size_t>> &reg, const std::vector<value_type> &alphas,
                    const std::vector<std::reference_wrapper<const AR>> &xx,
                    std::vector<std::reference_wrapper<AL>> &yy) override {
    return fused_axpyAny(reg, alphas, xx, yy);
  }
  void fused_axpyRR(const std::vector<std::tuple<size_t, size_t, size_t>> &reg, const std::vector<value_type> &alphas,
                    const std::vector<std::reference_wrapper<const AR>> &xx,
                    std::vector<std::reference_wrapper<AR>> &yy) override {
    return fused_axpyAny(reg, alphas, xx, yy);
  }

  template <typename T, typename S>
  void fused_dotAny(const std::vector<std::tuple<size_t, size_t, size_t>> &reg,
                    const std::vector<std::reference_wrapper<const T>> &xx,
                    const std::vector<std::reference_wrapper<const S>> &yy,
                    std::vector<std::reference_wrapper<value_type>> &out) {
    for (size_t i = 0; i < reg.size(); ++i) {
      size_t xi, yi, zi;
      std::tie(xi, yi, zi) = reg[i];
      out[i].get() = dot(xx[xi].get(), yy[yi].get());
    }
  }
  void fused_dotLL(const std::vector<std::tuple<size_t, size_t, size_t>> &reg,
                   const std::vector<std::reference_wrapper<const AL>> &xx,
                   const std::vector<std::reference_wrapper<const AL>> &yy,
                   std::vector<std::reference_wrapper<value_type>> &out) override {
    return fused_dotAny<AL, AL>(reg, xx, yy, out);
  }
  void fused_dotLR(const std::vector<std::tuple<size_t, size_t, size_t>> &reg,
                   const std::vector<std::reference_wrapper<const AL>> &xx,
                   const std::vector<std::reference_wrapper<const AR>> &yy,
                   std::vector<std::reference_wrapper<value_type>> &out) override {
    return fused_dotAny(reg, xx, yy, out);
  }
  void fused_dotRL(const std::vector<std::tuple<size_t, size_t, size_t>> &reg,
                   const std::vector<std::reference_wrapper<const AR>> &xx,
                   const std::vector<std::reference_wrapper<const AL>> &yy,
                   std::vector<std::reference_wrapper<value_type>> &out) override {
    return fused_dotAny(reg, xx, yy, out);
  }
  void fused_dotRR(const std::vector<std::tuple<size_t, size_t, size_t>> &reg,
                   const std::vector<std::reference_wrapper<const AR>> &xx,
                   const std::vector<std::reference_wrapper<const AR>> &yy,
                   std::vector<std::reference_wrapper<value_type>> &out) override {
    return fused_dotAny(reg, xx, yy, out);
  }
}; // namespace linalg

} // namespace array
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERITERABLE_H
