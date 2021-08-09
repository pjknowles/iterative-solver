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
#include <molpro/Profiler.h>
#include <future>
#include <molpro/cblas.h>

using molpro::linalg::itsolv::VecRef;
using molpro::linalg::itsolv::CVecRef;
using molpro::linalg::itsolv::subspace::Matrix;

namespace molpro::linalg::array::util {

template <class AL, class AR = AL>
void gemm_outer_distr_distr(const Matrix<typename array::mapped_or_value_type_t<AL>> alphas, const CVecRef<AR> &xx,
                            const VecRef<AL> &yy) {
  auto prof = molpro::Profiler::single();
  prof->push("gemm_outer_distr_distr (unbuffered)");
  *prof += alphas.rows() * alphas.cols() * yy[0].get().local_buffer()->size() * 2;
  for (size_t ii = 0; ii < alphas.rows(); ++ii) {
    auto loc_x = xx.at(ii).get().local_buffer();
    for (size_t jj = 0; jj < alphas.cols(); ++jj) {
      auto loc_y = yy[jj].get().local_buffer();
      for (size_t i = 0; i < loc_y->size(); ++i)
        (*loc_y)[i] += alphas(ii, jj) * (*loc_x)[i];
    }
  }
}
/**
template <class AL>
void gemm_outer_distr_distr(const Matrix<typename array::mapped_or_value_type_t<AL>> alphas,
                            const CVecRef<DistrArrayFile> &xx,
                            const VecRef<AL> &yy) {
  if (alphas.rows() > xx.size()){
    throw std::out_of_range("gemm_outer_distr_distr: dimensions of xx and alpha are different.");
  }
  std::cout << "the cool gemm_outer_outer\n";
  auto prof = molpro::Profiler::single();
  prof->start("gemm_outer_distr_distr (buffered)");
  *prof += alphas.rows() * alphas.cols() * yy[0].get().local_buffer()->size() * 2;

  std::vector<BufferManager> buffers;
  std::vector<BufferManager::Iterator> buffer_iterators;
  buffers.reserve(xx.size());
  buffer_iterators.reserve(xx.size());
  // TODO: reserve size of these vectors
  for (size_t j = 0; j < xx.size(); ++j){
    buffers.emplace_back(BufferManager(xx.at(j).get()));
    buffer_iterators.emplace_back(buffers[j].begin());
  }

  const int buf_stride = buffers[1].get_array_ptr() - buffers[0].get_array_ptr();
  const int yy_stride = yy[1].get().local_buffer()->data() - yy[0].get().local_buffer()->data();

  int previous_buf_stride;
  int previous_yy_stride;

  // check contiguity
  for (size_t j = 0; j<xx.size()-1; ++j){
    int curr_buf_stride = buffers[j+1].get_array_ptr() - buffers[j].get_array_ptr();
    int curr_yy_stride = yy[j+1].get().local_buffer()->data() - yy[j].get().local_buffer()->data();
    //std::cout << curr_buf_stride << ", " << curr_yy_stride << "\n";
    if (j<0){
      if (curr_yy_stride != previous_yy_stride){
        throw std::runtime_error("yy_stride is not contiguous\n");
      }
      if (curr_buf_stride != previous_buf_stride){
        throw std::runtime_error("buf_stride is not contiguous\n");
      }
    }
    previous_buf_stride = buf_stride;
    previous_yy_stride = yy_stride;
  } 

  size_t chunk_size = buffers[0].chunk_size;
  for (size_t curr_chunk = 0; curr_chunk < alphas.cols(); curr_chunk += chunk_size){

    auto buf_rows = buffer_iterators[0]->size();

    const int xx_cols = alphas.rows();
    const int alphas_rows = alphas.rows();
    const int alphas_cols = alphas.cols();
    const int yy_cols = yy.size();
    const int beta=1;
    const int alpha=1;

    const int lda = buf_stride;
    const int ldb = 1;
    const int ldc = yy_stride;

    // todo:
      // build a very detailed test
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, buf_rows, yy_cols, xx_cols, alpha, (*buffer_iterators[0]).data(),
                lda, alphas.data().data(), ldb, 1, yy[curr_chunk].get().local_buffer()->data(), ldc);

    for (size_t j = 1; j < xx.size(); ++j) { 
      ++buffer_iterators[j];
    }

  }
  prof->stop();
}
*/

template <class AL>
void gemm_outer_distr_distr(const Matrix<typename array::mapped_or_value_type_t<AL>> alphas,
                            const CVecRef<DistrArrayFile> &xx,
                            const VecRef<AL> &yy) {
  if (alphas.rows() > xx.size()){
    throw std::out_of_range("gemm_outer_distr_distr: dimensions of xx and alpha are different.");
  }
  auto prof = molpro::Profiler::single();
  prof->start("gemm_outer_distr_distr (buffered)");
  *prof += alphas.rows() * alphas.cols() * yy[0].get().local_buffer()->size() * 2;
  for (size_t ii = 0; ii < alphas.rows(); ++ii) {
    BufferManager x_buf = BufferManager(xx.at(ii).get());
    size_t offset = 0;
    for (auto buffer = x_buf.begin(); buffer != x_buf.end(); offset += x_buf.chunk_size, ++buffer) {
      size_t jj;
      for (jj = 0; jj < alphas.cols(); ++jj) { 
        auto loc_y = yy[jj].get().local_buffer();
        for (size_t i = 0; i < x_buf.chunk_size && i + offset < loc_y->size(); ++i){ 
          (*loc_y)[i + offset]  += alphas(ii, jj) * (*buffer)[i];
        }
      }
    }
  }
  prof->stop();
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
  //const size_t spacing = 1;
  using value_type = typename array::mapped_or_value_type_t<AL>;
  auto mat = Matrix<value_type>({xx.size(), yy.size()});
  if (xx.size() == 0 || yy.size() == 0) return mat;
  auto prof = molpro::Profiler::single()->push("gemm_inner_distr_distr (unbuffered)");
  prof += mat.cols() * mat.rows() * xx.at(0).get().local_buffer()->size()*2;
  for (size_t j = 0; j < mat.cols(); ++j) {
    auto loc_y = yy.at(j).get().local_buffer();
    for (size_t i = 0; i < mat.rows(); ++i) {
      auto loc_x = xx.at(i).get().local_buffer();
      mat(i, j) = std::inner_product(begin(*loc_x), end(*loc_x), begin(*loc_y), (value_type)0);
      //mat(i,j) += cblas_ddot(end(*loc_x) - begin(*loc_x), begin(*loc_x), spacing, begin(*loc_y), spacing);
    }
  }
#ifdef HAVE_MPI_H
  MPI_Allreduce(MPI_IN_PLACE, const_cast<value_type *>(mat.data().data()), mat.size(), MPI_DOUBLE, MPI_SUM,
                xx.at(0).get().communicator());
#endif
  return mat;
}

template <class AL> // async loading of loc_y into a managed buffer with buffers pre-instantiated
Matrix<typename array::mapped_or_value_type_t<AL>> gemm_inner_distr_distr(const CVecRef<AL> &xx,
                                                                          const CVecRef<DistrArrayFile> &yy) {
  using value_type = typename array::mapped_or_value_type_t<AL>;
  const size_t spacing = 1;
  auto mat = Matrix<value_type>({xx.size(), yy.size()});
  mat.fill(0);
  if (xx.size() == 0 || yy.size() == 0) return mat;
  auto prof = molpro::Profiler::single()->push("gemm_inner_distr_distr (buffered)");
  prof += mat.cols() * mat.rows() * xx.at(0).get().local_buffer()->size()*2;

  std::vector<BufferManager> buffers;
  for (size_t j = 0; j < mat.cols(); ++j){
    buffers.emplace_back(BufferManager(yy.at(j).get()));
  }

  for (size_t j = 0; j < mat.cols(); ++j) {
    size_t offset = 0;
    for (auto buffer = buffers[j].begin(); buffer != buffers[j].end(); offset += buffers[j].chunk_size, ++buffer) { 
      for (size_t i = 0; i < mat.rows(); ++i) {
        size_t buflen = buffer->cend() - buffer->cbegin();
        mat(i,j) += cblas_ddot(buflen, buffer->cbegin(), spacing, xx.at(i).get().local_buffer()->data() + offset,
                              spacing);
        //mat(i, j) += std::inner_product(buffer->cbegin(), buffer->cend(), xx.at(i).get().local_buffer()->data() + offset, (value_type)0);
      }
    }
  }
#ifdef HAVE_MPI_H
  MPI_Allreduce(MPI_IN_PLACE, const_cast<value_type *>(mat.data().data()), mat.size(), MPI_DOUBLE, MPI_SUM,
                xx.at(0).get().communicator());
#endif
  return mat;
}

/**
template <class AL>
Matrix<typename array::mapped_or_value_type_t<AL>> gemm_inner_distr_distr(const CVecRef<AL> &xx,
                                                                          const CVecRef<DistrArrayFile> &yy) {
  using value_type = typename array::mapped_or_value_type_t<AL>;
  auto mat = Matrix<value_type>({xx.size(), yy.size()});
  if (xx.size() == 0 || yy.size() == 0) return mat;
  auto prof = molpro::Profiler::single()->push("gemm_inner_distr_distr (buffered)");
  prof += mat.cols() * mat.rows() * xx.at(0).get().local_buffer()->size()*2;
  auto loc_y = yy.at(0).get().local_buffer();
  auto loc_y_future = std::async(std::launch::async, [yy] { return yy.at(1).get().local_buffer(); });
  for (size_t j = 0; j < mat.cols(); ++j) {
    if (j>0){
      loc_y = loc_y_future.get();
    }
    if (j < mat.cols()-1){
      loc_y_future = std::async(std::launch::async, [yy,j] { return yy.at(j+1).get().local_buffer(); });
    }
    for (size_t i = 0; i < mat.rows(); ++i) {
      auto loc_x = xx.at(i).get().local_buffer();
      mat(i, j) = std::inner_product(begin(*loc_x), end(*loc_x), begin(*loc_y), (value_type)0);
      //mat(i,j) += cblas_ddot(loc_x->size(), begin(*loc_x), spacing, begin(*loc_y), spacing);
    }
  }
#ifdef HAVE_MPI_H
  MPI_Allreduce(MPI_IN_PLACE, const_cast<value_type *>(mat.data().data()), mat.size(), MPI_DOUBLE, MPI_SUM,
                xx.at(0).get().communicator());
#endif
  return mat;
}
*/
/**
template <class AL> asyncronous loading of loc_y into a buffer managed by buffermanager
Matrix<typename array::mapped_or_value_type_t<AL>> gemm_inner_distr_distr(const CVecRef<AL> &xx,
                                                                          const CVecRef<DistrArrayFile> &yy) {
  using value_type = typename array::mapped_or_value_type_t<AL>;
  auto mat = Matrix<value_type>({xx.size(), yy.size()});
  mat.fill(0);
  if (xx.size() == 0 || yy.size() == 0) return mat;
  auto prof = molpro::Profiler::single()->push("gemm_inner_distr_distr (buffered)");
  prof += mat.cols() * mat.rows() * xx.at(0).get().local_buffer()->size()*2;
  for (size_t j = 0; j < mat.cols(); ++j) {
    BufferManager y_buf = BufferManager(yy.at(j).get());
    size_t offset = 0;
    for (auto buffer = y_buf.begin(); buffer != y_buf.end(); offset += y_buf.chunk_size, ++buffer) { 
      for (size_t i = 0; i < mat.rows(); ++i) {
        mat(i, j) += std::inner_product(buffer->cbegin(), buffer->cend(), xx.at(i).get().local_buffer()->data() + offset, (value_type)0);
      }
    }
  }
#ifdef HAVE_MPI_H
  MPI_Allreduce(MPI_IN_PLACE, const_cast<value_type *>(mat.data().data()), mat.size(), MPI_DOUBLE, MPI_SUM,
                xx.at(0).get().communicator());
#endif
  return mat;
}
*/
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
