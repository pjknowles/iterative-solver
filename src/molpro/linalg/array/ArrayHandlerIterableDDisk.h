#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERITERABLEDDISK_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERITERABLEDDISK_H
#include <algorithm>
#include <iterator>
#include <molpro/linalg/array/ArrayHandler.h>
#include <molpro/linalg/array/ArrayHandlerIterable.h>
#include <molpro/linalg/array/Span.h>
#include <molpro/linalg/array/util/gemm.h>
#include <molpro/linalg/array/util/select.h>
#include <molpro/linalg/scalar_traits.h>

namespace molpro::linalg::array {

namespace util {
//! A view of all the elements of a contiguous array
template <class R>
auto contiguous_span(R &r) {
  return Span<std::remove_pointer_t<decltype(std::data(r))>>(std::data(r), std::size(r));
}
} // namespace util

/*!
 * @brief Array handler between a contiguous array (left), e.g. std::vector or Span, and a distributed array (right),
 * typically one held on disk.
 *
 * The contiguous array holds all elements on every process, i.e. it is replicated, while each process works on its
 * local section of the distributed array, reading it page by page through get() if it is on disk. The contiguous
 * array is addressed in place. Operations that write into the contiguous array finish by replicating it over all
 * processes. When the distributed array lives on a single process (e.g. molpro::mpi::comm_self()), there is no
 * communication at all.
 *
 * ArrayHandlerDDiskIterable handles the same pair of arrays the other way round.
 */
template <class AL, class AR>
class ArrayHandlerIterableDDisk : public ArrayHandler<AL, AR> {
public:
  using typename ArrayHandler<AL, AR>::value_type_L;
  using typename ArrayHandler<AL, AR>::value_type_R;
  using typename ArrayHandler<AL, AR>::value_type;
  using typename ArrayHandler<AL, AR>::value_type_abs;
  using typename ArrayHandler<AL, AR>::ProxyHandle;
  static_assert(is_contiguous_v<AL>, "ArrayHandlerIterableDDisk: the left array must be contiguous");
  static_assert(std::is_same_v<value_type, typename AR::value_type>,
                "ArrayHandlerIterableDDisk: the arrays must have the same value type");

  ProxyHandle lazy_handle() override { return this->lazy_handle(*this); };

  using ArrayHandler<AL, AR>::lazy_handle;
  using ArrayHandler<AL, AR>::error;

  AL copy(const AR &source) override {
    auto result = util::allocate_like<AL>(source.size());
    copy(result, source);
    return result;
  };

  void copy(AL &x, const AR &y) override {
    this->m_counter->copy++;
    y.copy_to(util::contiguous_span(x));
  };

  void scal(value_type alpha, AL &x) override {
    this->m_counter->scal++;
    for (auto &el : x)
      el *= alpha;
  }

  void fill(value_type alpha, AL &x) override { std::fill(std::begin(x), std::end(x), alpha); }

  void axpy(value_type alpha, const AR &x, AL &y) override {
    this->m_counter->axpy++;
    x.axpy_to(alpha, util::contiguous_span(y));
  }

  //! The hermitian inner product <x|y>, i.e. the conjugate of <y|x>
  value_type dot(const AL &x, const AR &y) override {
    this->m_counter->dot++;
    return molpro::linalg::conjugate(y.dot(util::contiguous_span(x)));
  }

  void gemm_outer(const Matrix<value_type> alphas, const CVecRef<AR> &xx, const VecRef<AL> &yy) override {
    this->m_counter->gemm_outer++;
    util::gemm_outer_distr_contiguous(alphas, xx, yy);
  }

  Matrix<value_type> gemm_inner(const CVecRef<AL> &xx, const CVecRef<AR> &yy) override {
    this->m_counter->gemm_inner++;
    return util::gemm_inner_contiguous_distr(xx, yy);
  }

  std::map<size_t, value_type_abs> select_max_dot(size_t n, const AL &x, const AR &y) override {
    return y.select_max_dot(n, util::contiguous_span(x));
  }

  std::map<size_t, value_type> select(size_t n, const AL &x, bool max = false, bool ignore_sign = false) override {
    if (n > std::size(x))
      error("ArrayHandlerIterableDDisk::select() n is too large");
    return util::select<AL, value_type>(n, x, max, ignore_sign);
  }
};

} // namespace molpro::linalg::array

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERITERABLEDDISK_H
