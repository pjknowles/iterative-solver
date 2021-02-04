#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERSPARSE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERSPARSE_H
#include <molpro/linalg/array/ArrayHandler.h>
#include <molpro/linalg/array/util/select_max_dot.h>
#include <molpro/linalg/array/util/gemm.h>

using molpro::linalg::array::util::gemm_outer_default;
using molpro::linalg::array::util::gemm_inner_default;

namespace molpro::linalg::array {

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

  void copy(AL &x, const AR &y) override {
    using std::begin;
    using std::end;
    x.clear();
    std::copy(begin(y), end(y), std::inserter(x, x.begin()));
  };

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
  
  void gemm_outer(const Matrix<value_type> alphas, const CVecRef<AR> &xx, const VecRef<AL> &yy) override {
    gemm_outer_default(*this, alphas, xx, yy);
  }
  
  Matrix<value_type> gemm_inner(const CVecRef<AL> &xx, const CVecRef<AR> &yy) override {
    return gemm_inner_default(*this, xx, yy);
  }
  
  std::map<size_t, value_type_abs> select_max_dot(size_t n, const AL &x, const AR &y) override {
    if (n > x.size() || n > y.size())
      error("ArrayHandlerSparse::select_max_dot() n is too large");
    return util::select_max_dot_sparse<AL, AR, value_type, value_type_abs>(n, x, y);
  }

  ProxyHandle lazy_handle() override { return this->lazy_handle(*this); };

protected:
  using ArrayHandler<AL, AR>::error;
  using ArrayHandler<AL, AR>::lazy_handle;
};

} // namespace molpro::linalg::array

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERSPARSE_H
