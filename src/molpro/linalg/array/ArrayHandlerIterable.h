#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERITERABLE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERITERABLE_H
#include "molpro/linalg/array/ArrayHandler.h"
#include <numeric>

namespace molpro {
namespace linalg {
namespace array {

//! Array handler for two containers both of which can be iterated through using begin()/cbegin() and end()/cend()
//! member functions.
template <typename AL, typename AR = AL> class ArrayHandlerIterable : public ArrayHandler<AL, AR> {
public:
  using typename ArrayHandler<AL, AR>::value_type;
  ArrayHandlerIterable() = default;

  value_type dot(const AL &x, const AL &y) override {
    if (x.size() > y.size())
      error("ArrayHandlerIterable::dot() incompatible x and y arrays, x.size() > y.size()");
    return std::inner_product(x.cbegin(), x.cend(), y.cbegin(), 0);
  }

  AL &axpy(AL &x, value_type alpha, const AR &y) override {
    if (x.size() > y.size())
      error("ArrayHandlerIterable::dot() incompatible x and y arrays, x.size() > y.size()");
    auto iy = y.cbegin();
    for (auto ix = x.begin(); ix != x.end(); ++ix, ++iy)
      *ix += alpha * (*iy);
    return x;
  }

  void error(std::string message) override { throw std::runtime_error(message); }

protected:
  template <typename T, typename S> T copyAny(const S &source) {
    T result{};
    std::copy(source.cbegin(), source.cend(), std::back_inserter(result));
    return result;
  }
  AL copyL(const AR &source) override { return copyAny<AL, AR>(source); }
  AR copyR(const AR &source) override { return copyAny<AR, AL>(source); }

  void fused_axpy(const std::vector<std::pair<size_t, size_t>> &reg, std::vector<std::reference_wrapper<AL>> &xx,
                  std::vector<std::reference_wrapper<const AR>> &yy, std::vector<value_type> alphas) override {
    for (size_t i = 0; i < reg.size(); ++i) {
      size_t xi, yi;
      std::tie(xi, yi) = reg[i];
      axpy(xx[xi].get(), alphas[i], yy[yi].get());
    }
  } // namespace array

  void fused_dot(const std::vector<std::pair<size_t, size_t>> &reg, std::vector<std::reference_wrapper<const AL>> &xx,
                 std::vector<std::reference_wrapper<const AR>> &yy,
                 std::vector<std::reference_wrapper<value_type>> &out) override {
    for (size_t i = 0; i < reg.size(); ++i) {
      size_t xi, yi;
      std::tie(xi, yi) = reg[i];
      out[i].get() = dot(xx[xi].get(), yy[yi].get());
    }
  }
}; // namespace linalg

} // namespace array
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERITERABLE_H
