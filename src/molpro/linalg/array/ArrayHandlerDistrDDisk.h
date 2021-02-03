#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERDISTRDDISK_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERDISTRDDISK_H

#include <map>
#include <molpro/linalg/array/ArrayHandler.h>

namespace molpro::linalg::array {

template <class AL, class AR = AL>
class ArrayHandlerDistrDDisk : public ArrayHandler<AL, AR> {
public:
  using typename ArrayHandler<AL, AR>::value_type_L;
  using typename ArrayHandler<AL, AR>::value_type_R;
  using typename ArrayHandler<AL, AR>::value_type;
  using typename ArrayHandler<AL, AR>::value_type_abs;
  using typename ArrayHandler<AL, AR>::ProxyHandle;

  ProxyHandle lazy_handle() override { return this->lazy_handle(*this); };

  using ArrayHandler<AL, AR>::lazy_handle;
  using ArrayHandler<AL, AR>::error;

  AL copy(const AR &source) override { return AL{source}; };

  void copy(AL &x, const AR &y) override { x.copy(y); };

  void scal(value_type alpha, AL &x) override { x.scal(alpha); }

  void fill(value_type alpha, AL &x) override { x.fill(alpha); }

  void axpy(value_type alpha, const AR &x, AL &y) override { y.axpy(alpha, x); }

  value_type dot(const AL &x, const AR &y) override { return x.dot(y); }
  
  void gemm_outer(const Matrix<value_type> alphas, const CVecRef<AR> &xx, const VecRef<AL> &yy) override {
    for (size_t ii = 0; ii < alphas.cols(); ++ii) {
      auto loc_x = xx.at(ii).get().local_buffer();
      for (size_t jj = 0; jj < alphas.rows(); ++jj) {
        auto loc_y = yy[jj].get().local_buffer();
        for (size_t i = 0; i < loc_x->size(); ++i)
          (*loc_y)[i] += alphas(ii, jj) * (*loc_x)[i];
      }
    }
  }
  
  Matrix<value_type> gemm_inner(const CVecRef<AL> &xx, const CVecRef<AR> &yy) override {
    auto mat = Matrix<value_type>({xx.size(), yy.size()});
    for (size_t i = 0; i < mat.cols(); ++i) {
      auto loc_y = yy.at(i).get().local_buffer();
      for (size_t j = 0; j < mat.rows(); ++j) {
        auto loc_x = xx.at(j).get().local_buffer();
        mat(i, j) = std::inner_product(begin(*loc_x), end(*loc_x), begin(*loc_y), (value_type)0);
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, const_cast<value_type *>(mat.data().data()), mat.size(), MPI_DOUBLE, MPI_SUM,
                  xx[0].get().communicator());
    return mat;
  }
  
  std::map<size_t, value_type_abs> select_max_dot(size_t n, const AL &x, const AR &y) override {
    return x.select_max_dot(n, y);
  }
};

} // namespace molpro::linalg::array

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERDISTRDDISK_H
