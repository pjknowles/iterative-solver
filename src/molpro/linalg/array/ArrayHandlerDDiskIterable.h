#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERDDISKITERABLE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERDDISKITERABLE_H
#include <molpro/linalg/array/ArrayHandlerIterableDDisk.h>
#include <molpro/linalg/array/ArrayHandlerWithCopyFunc.h>
#include <molpro/linalg/array/util/gemm.h>
#include <molpro/mpi.h>
#include <string>
#include <type_traits>

namespace molpro::linalg::array {

/*!
 * @brief Array handler between a distributed array (left), typically one held on disk, and a contiguous array
 * (right), e.g. std::vector or Span.
 *
 * This is ArrayHandlerIterableDDisk the other way round, see there for how the arrays are accessed.
 *
 * A new distributed array made as a copy of a contiguous one is created by a copy function. By default, the new array
 * is constructed from its size, molpro::mpi::comm_self() (the contiguous array is not distributed) and a directory.
 */
template <class AL, class AR>
class ArrayHandlerDDiskIterable : public ArrayHandlerWithCopyFunc<AL, AR> {
public:
  using typename ArrayHandler<AL, AR>::value_type_L;
  using typename ArrayHandler<AL, AR>::value_type_R;
  using typename ArrayHandler<AL, AR>::value_type;
  using typename ArrayHandler<AL, AR>::value_type_abs;
  using typename ArrayHandler<AL, AR>::ProxyHandle;
  static_assert(is_contiguous_v<AR>, "ArrayHandlerDDiskIterable: the right array must be contiguous");
  static_assert(std::is_same_v<value_type, typename AL::value_type>,
                "ArrayHandlerDDiskIterable: the arrays must have the same value type");

  //! New arrays are created on molpro::mpi::comm_self(), with any file in directory
  explicit ArrayHandlerDDiskIterable(const std::string &directory = ".")
      : ArrayHandlerWithCopyFunc<AL, AR>([directory](const AR &source) -> AL {
          if constexpr (std::is_constructible_v<AL, size_t, MPI_Comm, std::string>) {
            auto result = AL(std::size(source), molpro::mpi::comm_self(), directory);
            result.copy(util::contiguous_span(source));
            return result;
          } else
            throw std::logic_error("ArrayHandlerDDiskIterable: cannot construct the distributed array from its size, "
                                   "communicator and directory; a copy function must be given");
        }) {}

  //! Constructor taking a copy function, which returns a new distributed array as a copy of a contiguous one
  explicit ArrayHandlerDDiskIterable(std::function<AL(const AR &)> copy_func)
      : ArrayHandlerWithCopyFunc<AL, AR>(std::move(copy_func)) {}

  ProxyHandle lazy_handle() override { return this->lazy_handle(*this); };

  using ArrayHandler<AL, AR>::lazy_handle;
  using ArrayHandler<AL, AR>::error;
  using ArrayHandlerWithCopyFunc<AL, AR>::copy;

  void copy(AL &x, const AR &y) override {
    this->m_counter->copy++;
    x.copy(util::contiguous_span(y));
  };

  void scal(value_type alpha, AL &x) override {
    this->m_counter->scal++;
    x.scal(alpha);
  }

  void fill(value_type alpha, AL &x) override { x.fill(alpha); }

  void axpy(value_type alpha, const AR &x, AL &y) override {
    this->m_counter->axpy++;
    y.axpy(alpha, util::contiguous_span(x));
  }

  value_type dot(const AL &x, const AR &y) override {
    this->m_counter->dot++;
    return x.dot(util::contiguous_span(y));
  }

  void gemm_outer(const Matrix<value_type> alphas, const CVecRef<AR> &xx, const VecRef<AL> &yy) override {
    this->m_counter->gemm_outer++;
    util::gemm_outer_contiguous_distr(alphas, xx, yy);
  }

  Matrix<value_type> gemm_inner(const CVecRef<AL> &xx, const CVecRef<AR> &yy) override {
    this->m_counter->gemm_inner++;
    return util::gemm_inner_distr_contiguous(xx, yy);
  }

  std::map<size_t, value_type_abs> select_max_dot(size_t n, const AL &x, const AR &y) override {
    return x.select_max_dot(n, util::contiguous_span(y));
  }

  std::map<size_t, value_type> select(size_t n, const AL &x, bool max = false, bool ignore_sign = false) override {
    return x.select(n, max, ignore_sign);
  }
};

} // namespace molpro::linalg::array

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERDDISKITERABLE_H
