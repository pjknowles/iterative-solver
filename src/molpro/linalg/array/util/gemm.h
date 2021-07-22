#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_GEMM_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_GEMM_H
#ifdef HAVE_MPI_H
#include <mpi.h>
#endif
#include <iostream>
#include <vector>
#include <numeric>
#include <molpro/linalg/array/type_traits.h>
#include <molpro/linalg/itsolv/wrap_util.h>
#include <molpro/linalg/itsolv/subspace/Matrix.h>
#include <molpro/linalg/array/DistrArrayFile.h>

using molpro::linalg::itsolv::VecRef;
using molpro::linalg::itsolv::CVecRef;
using molpro::linalg::itsolv::subspace::Matrix;

namespace molpro::linalg::array::util {

template <class AL, class AR = AL>
void gemm_outer_distr_distr(const Matrix<typename array::mapped_or_value_type_t<AL>> alphas, const CVecRef<AR> &xx,
                            const VecRef<AL> &yy) {
  for (size_t ii = 0; ii < alphas.rows(); ++ii) {
    auto loc_x = xx.at(ii).get().local_buffer();
    for (size_t jj = 0; jj < alphas.cols(); ++jj) {
      auto loc_y = yy[jj].get().local_buffer();
      for (size_t i = 0; i < loc_y->size(); ++i)
        (*loc_y)[i] += alphas(ii, jj) * (*loc_x)[i];
    }
  }
}

template <class AL>
void gemm_outer_distr_distr(const Matrix<typename array::mapped_or_value_type_t<AL>> alphas,
                            const CVecRef<DistrArrayFile> &xx,
                            const VecRef<AL> &yy) {
  for (size_t ii = 0; ii < alphas.rows(); ++ii) { // should be more outermost
    BufferManager x_buf = BufferManager(xx.at(ii).get());
    size_t offset = 0;
    for (auto buffer = x_buf.begin(); buffer != x_buf.end(); offset += x_buf.chunk_size, ++buffer) { // should be outermost
      size_t jj;
      for (jj = 0; jj < alphas.cols(); ++jj) { 
        auto loc_y = yy[jj].get().local_buffer();
        for (size_t i = 0; i < loc_y->size(); ++i){ // should be middle
          (*loc_y)[i + offset]  += alphas(ii, jj) * (*buffer)[i];
        }
      }
    }
  }
}

template <class AL, class AR = AL>
void gemm_outer_distr_sparse(const Matrix<typename array::mapped_or_value_type_t<AL>> alphas, const CVecRef<AR> &xx,
                             const VecRef<AL> &yy) {
  for (size_t ii = 0; ii < alphas.cols(); ++ii) {
    auto loc_y = yy[ii].get().local_buffer();
    for (size_t jj = 0; jj < alphas.rows(); ++jj) {
      if (loc_y->size() > 0) {
        size_t i;
        typename array::mapped_or_value_type_t<AL> v;
        for (auto it = xx.at(jj).get().lower_bound(loc_y->start());
             it != xx.at(jj).get().upper_bound(loc_y->start() + loc_y->size() - 1); ++it) {
          std::tie(i, v) = *it;
          (*loc_y)[i - loc_y->start()] += alphas(jj, ii) * v;
        }
      }
    }
  }
}

template <class Handler, class AL, class AR = AL>
void gemm_outer_default(Handler &handler, const Matrix<typename Handler::value_type> alphas,
                        const CVecRef<AR> &xx, const VecRef<AL> &yy) {
  for (size_t ii = 0; ii < alphas.rows(); ++ii) {
    for (size_t jj = 0; jj < alphas.cols(); ++jj) {
      handler.axpy(alphas(ii, jj), xx.at(ii).get(), yy[jj].get());
    }
  }
}

template <class AL, class AR = AL>
Matrix<typename array::mapped_or_value_type_t<AL>> gemm_inner_distr_distr(const CVecRef<AL> &xx,
                                                                          const CVecRef<AR> &yy) {
  using value_type = typename array::mapped_or_value_type_t<AL>;
  auto mat = Matrix<value_type>({xx.size(), yy.size()});
  if (xx.size() == 0 || yy.size() == 0) return mat;
  for (size_t j = 0; j < mat.cols(); ++j) {
    auto loc_y = yy.at(j).get().local_buffer();
    for (size_t i = 0; i < mat.rows(); ++i) {
      auto loc_x = xx.at(i).get().local_buffer();
      mat(i, j) = std::inner_product(begin(*loc_x), end(*loc_x), begin(*loc_y), (value_type)0);
    }
  }
#ifdef HAVE_MPI_H
  MPI_Allreduce(MPI_IN_PLACE, const_cast<value_type *>(mat.data().data()), mat.size(), MPI_DOUBLE, MPI_SUM,
                xx.at(0).get().communicator());
#endif
  return mat;
}

template <class AL, class std::enable_if<!std::is_same<AL,DistrArrayFile>::value>::type* = nullptr>
Matrix<typename array::mapped_or_value_type_t<AL>> gemm_inner_distr_distr(const CVecRef<AL> &xx,
                                                                          const CVecRef<DistrArrayFile> &yy) {
  if (xx == yy){
    throw std::invalid_argument("Cannot gemm a VecRef with itself.");
  }
  using value_type = typename array::mapped_or_value_type_t<AL>;
  auto mat = Matrix<value_type>({xx.size(), yy.size()});
  mat.fill(0);
  if (xx.size() == 0 || yy.size() == 0) return mat;
  for (size_t j = 0; j < mat.cols(); ++j) {
    BufferManager y_buf = BufferManager(yy.at(j).get());
    size_t offset = 0;
    for (auto buffer = y_buf.begin(); buffer != y_buf.end(); offset += y_buf.chunk_size, ++buffer) { // this needs to outside such that you don't read the file for every j
      for (size_t i = 0; i < mat.rows(); ++i) {
        mat(i, j) = std::inner_product(buffer->cbegin(), buffer->cend(), xx.at(i).get().local_buffer()->data() + offset, mat(i, j));
      }
    }
  }
#ifdef HAVE_MPI_H
  MPI_Allreduce(MPI_IN_PLACE, const_cast<value_type *>(mat.data().data()), mat.size(), MPI_DOUBLE, MPI_SUM,
                xx.at(0).get().communicator());
#endif
  return mat;
}

template <class AL, class AR = AL>
Matrix<typename array::mapped_or_value_type_t<AL>> gemm_inner_distr_sparse(const CVecRef<AL> &xx,
                                                                           const CVecRef<AR> &yy) {
  using value_type = typename array::mapped_or_value_type_t<AL>;
  auto mat = Matrix<value_type>({xx.size(), yy.size()});
  if (xx.size() == 0 || yy.size() == 0) return mat;
  for (size_t i = 0; i < mat.rows(); ++i) {
    auto loc_x = xx.at(i).get().local_buffer();
    for (size_t j = 0; j < mat.cols(); ++j) {
      mat(i, j) = 0;
      if (loc_x->size() > 0) {
        size_t k;
        value_type v;
        for (auto it = yy.at(j).get().lower_bound(loc_x->start());
             it != yy.at(j).get().upper_bound(loc_x->start() + loc_x->size() - 1); ++it) {
          std::tie(k, v) = *it;
          mat(i, j) += (*loc_x)[k - loc_x->start()] * v;
        }
      }
    }
  }
#ifdef HAVE_MPI_H
  MPI_Allreduce(MPI_IN_PLACE, const_cast<value_type*>(mat.data().data()), mat.size(), MPI_DOUBLE, MPI_SUM,
                xx.at(0).get().communicator());
#endif
  return mat;
}

template <class Handler, class AL, class AR = AL>
Matrix<typename Handler::value_type> gemm_inner_default(Handler &handler, const CVecRef<AL> &xx, const CVecRef<AR> &yy) {
  auto mat = Matrix<typename Handler::value_type>({xx.size(), yy.size()});
  if (xx.size() == 0 || yy.size() == 0) return mat;
  for (size_t ii = 0; ii < mat.rows(); ++ii) {
    for (size_t jj = 0; jj < mat.cols(); ++jj) {
      mat(ii, jj) = handler.dot(xx.at(ii).get(), yy.at(jj).get());
    }
  }
  return mat;
}

} // namespace molpro::linalg::array::util

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_GEMM_H
