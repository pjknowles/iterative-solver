#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERDDISKSPARSE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERDDISKSPARSE_H
#include <map>
#include <molpro/linalg/array/ArrayHandler.h>
#include <numeric>

namespace molpro::linalg::array {
template <typename AL, typename AR, bool = has_mapped_type_v<AR>>
class ArrayHandlerDDiskSparse : public ArrayHandler<AL, AR> {};

/*!
 * @brief Array handler between an iterable and a sparse arrays. Iterable container must implement access operator.
 */
template <typename AL, typename AR>
class ArrayHandlerDDiskSparse<AL, AR, true> : public ArrayHandler<AL, AR> {
public:
  using typename ArrayHandler<AL, AR>::value_type_L;
  using typename ArrayHandler<AL, AR>::value_type_R;
  using typename ArrayHandler<AL, AR>::value_type;
  using typename ArrayHandler<AL, AR>::value_type_abs;
  using typename ArrayHandler<AL, AR>::ProxyHandle;

  AL copy(const AR &source) override {
    static_assert(true, "General copy from sparse to dense is ill-defined");
    return AL{};
  };
  void copy(AL &x, const AR &y) override { static_assert(true, "General copy from sparse to dense is ill-defined"); };

  void scal(value_type alpha, AL &x) override { static_assert(true, "Use ArrayHandlerDistr for unary operations"); };

  void fill(value_type alpha, AL &x) override { static_assert(true, "Use ArrayHandlerDistr for unary operations"); };

  void axpy(value_type alpha, const AR &x, AL &y) override { y.axpy(alpha, x); };

  value_type dot(const AL &x, const AR &y) override { return x.dot(y); };
  
  void gemm_outer(const Matrix<value_type> alphas, const CVecRef<AR> &xx, const VecRef<AL> &yy) override {
    for (size_t ii = 0; ii < alphas.cols(); ++ii) {
      auto loc_y = yy.at(ii).get().local_buffer();
      for (size_t jj = 0; jj < alphas.rows(); ++jj) {
        if (loc_y->size() > 0) {
          size_t i;
          value_type v;
          for (auto it = xx.at(jj).get().lower_bound(loc_y->start());
               it != xx.at(jj).get().upper_bound(loc_y->start() + loc_y->size() - 1); ++it) {
            std::tie(i, v) = *it;
            (*loc_y)[i - loc_y->start()] += alphas(ii, jj) * v;
          }
        }
      }
    }
  }
  
  Matrix<value_type> gemm_inner(const CVecRef<AL> &xx, const CVecRef<AR> &yy) override {
    auto mat = Matrix<value_type>({xx.size(), yy.size()});
    for (size_t i = 0; i < mat.cols(); ++i) {
      auto loc_x = xx.at(i).get().local_buffer();
      for (size_t j = 0; j < mat.rows(); ++j) {
        mat(i, j) = 0;
        if (loc_x->size() > 0) {
          size_t i;
          value_type v;
          for (auto it = yy.at(j).get().lower_bound(loc_x->start());
               it != yy.at(j).get().upper_bound(loc_x->start() + loc_x->size() - 1); ++it) {
            std::tie(i, v) = *it;
            mat(i, j) += (*loc_x)[i - loc_x->start()] * v;
          }
        }
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, const_cast<value_type*>(mat.data().data()), mat.size(), MPI_DOUBLE, MPI_SUM,
                  xx.at(0).get().communicator());
    return mat;
  }
  
  std::map<size_t, value_type_abs> select_max_dot(size_t n, const AL &x, const AR &y) override {
    return x.select_max_dot(n, y);
  }
  
  ProxyHandle lazy_handle() override { return this->lazy_handle(*this); };
  
protected:
  using ArrayHandler<AL, AR>::error;
  using ArrayHandler<AL, AR>::lazy_handle;
  using ArrayHandler<AL, AR>::m_lazy_handles;
};

} // namespace molpro::linalg::array

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERDDISKSPARSE_H
