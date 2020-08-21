#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERSPARSE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERSPARSE_H
#include <molpro/linalg/array/ArrayHandler.h>

namespace molpro {
namespace linalg {
namespace array {

/*!
 * @brief Array handler between two sparse arrays (e.g. std::map)
 */
template <typename AL, typename AR = AL>
class ArrayHandlerSparse : public ArrayHandler<AL, AR> {
public:
  using typename ArrayHandler<AL, AR>::value_type_L;
  using typename ArrayHandler<AL, AR>::value_type_R;
  using typename ArrayHandler<AL, AR>::value_type;
  using typename ArrayHandler<AL, AR>::value_type_abs;
  using typename ArrayHandler<AL, AR>::ProxyHandle;

  AL copy(const AR &source) override { return AL{source.begin(), source.end()}; };

  void scal(value_type alpha, AL &x) override {
    for (auto &el : x)
      el.second *= alpha;
  };

  void fill(value_type alpha, AL &x) override {
    for (auto &el : x)
      el.second = alpha;
  };

  void axpy(value_type alpha, const AR &x, AL &y) override {
    for (const auto &ix : x) {
      auto iy = y.find(ix.first);
      if (iy != y.end())
        iy->second += alpha * ix.second;
    }
  };

  value_type dot(const AL &x, const AR &y) override {
    auto tot = value_type{};
    for (const auto &ix : x) {
      const auto iy = y.find(ix.first);
      if (iy != y.end())
        tot += iy->second * ix.second;
    }
    return tot;
  };

  std::map<size_t, value_type_abs> select_max_dot(size_t n, const AL &x, const AR &y) override { return {}; }

  ProxyHandle lazy_handle() override { return this->lazy_handle(*this); };

protected:
  using ArrayHandler<AL, AR>::lazy_handle;
};

} // namespace array
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERSPARSE_H
