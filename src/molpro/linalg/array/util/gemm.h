#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_GEMM_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_GEMM_H
#ifdef HAVE_MPI_H
#include <mpi.h>
#endif
#include <future>
#include <iostream>
#include <molpro/Profiler.h>
#include <molpro/cblas.h>
#include <molpro/linalg/array/DistrArrayFile.h>
#include <molpro/linalg/array/type_traits.h>
#include <molpro/linalg/itsolv/subspace/Matrix.h>
#include <molpro/linalg/itsolv/wrap_util.h>
#include <numeric>
#include <vector>

using molpro::linalg::itsolv::CVecRef;
using molpro::linalg::itsolv::VecRef;
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

template <class AL>
void gemm_outer_distr_distr(const Matrix<typename array::mapped_or_value_type_t<AL>> alphas,
                            const CVecRef<DistrArrayFile> &xx, const VecRef<AL> &yy) {

  if (alphas.rows() != xx.size() ||
      alphas.cols() !=
          yy.size()) { // this should assert specifically that alphas.rows is xx.size AND alphas.cols is yy.size
    throw std::out_of_range("gemm_outer_distr_distr: dimensions of xx and alpha are different.");
  }

  if (alphas.rows() == 0 || alphas.cols() == 0) {
    return;
  }
  //
  //  if (alphas.rows() == 1 || alphas.cols() == 1){
  //    auto loc_x = xx.at(0).get().local_buffer();
  //    auto loc_y = yy[0].get().local_buffer();
  //    for (size_t i = 0; i < loc_y->size(); ++i)
  //      (*loc_y)[i] += alphas(0, 0) * (*loc_x)[i];
  //    return;
  //  }

  //  std::cout << "gemm_outer_distr_distr (buffered)\n";
  auto prof = molpro::Profiler::single();
  prof->start("gemm_outer_distr_distr (buffered)");
  *prof += alphas.rows() * alphas.cols() * yy[0].get().local_buffer()->size() * 2;

  bool yy_constant_stride = true;
  int previous_stride = 0;
  int yy_stride = yy.front().get().local_buffer()->size();
  for (size_t j = 0; j < std::max((size_t)1, yy.size()) - 1; ++j) {
    auto unique_ptr_j = yy.at(j).get().local_buffer()->data();
    auto unique_ptr_jp1 = yy.at(j + 1).get().local_buffer()->data();
    yy_stride = unique_ptr_jp1 - unique_ptr_j;
    if (j > 0)
      yy_constant_stride = yy_constant_stride && (yy_stride == previous_stride);
    previous_stride = yy_stride;
  }
  //if (not yy_constant_stride)
  //  throw std::runtime_error(
  //      "yy doesn't have a consistent stride, and general case not yet coded\n"); // TODO code up the non-constant
                                                                                  // stride case using dgemv

  const int buf_size =
      10; // The amount of memory allocated for buffering of xx. IN the case of double buffering, this means that the
          // actual buffers will be half this value. TODO: should be 8192 or similar large.
  // Eventual TODO: use a vector for this
  // std::vector<DistrArray::value_type> buffers_vec(buf_size*alphas.rows());
  // DistrArray::value_type* buffers_memory = buffers_vec.data();
  DistrArray::value_type *buffers_memory = new DistrArray::value_type[buf_size * alphas.rows()];
  std::vector<BufferManager> buffers;
  std::vector<BufferManager::Iterator> buffer_iterators;
  buffers.reserve(xx.size());
  buffer_iterators.reserve(xx.size());
  for (size_t j = 0; j < alphas.rows(); ++j) {
    buffers.emplace_back(
        BufferManager(xx.at(j).get(), buffers_memory + (buf_size * j), buf_size, BufferManager::buffertype::Double));
    buffer_iterators.emplace_back(buffers[j].begin());
  }
  const int N = alphas.cols(); // cols of alphas, cols of yy
  const int K = alphas.rows(); // cols of xx, rows of alphas
  size_t M;                    // = buffer_iterators.front()->size();
  const auto yy_data = yy[0].get().local_buffer()->data();
  const auto yy_size = yy[0].get().local_buffer()->size();
  for (size_t container_offset = 0; container_offset < yy_size; container_offset += M) {
    M = buffer_iterators.front()->size();
    if (false) {

      std::cout << "container loop start container_offset=" << container_offset << ", M=" << M
                << ", iterator size=" << buffer_iterators.front()->size() << std::endl;
      for (size_t i = 0; i < yy.size(); ++i) {
        std::cout << "presenting yy[" << i << "]";
        for (size_t j = 0; j < M; ++j)
          std::cout << " " << yy_data[container_offset + i * yy_stride + j];
        std::cout << std::endl;
      }
      for (size_t i = 0; i < xx.size(); ++i) {
        std::cout << "presenting xx[" << i << "]";
        for (size_t j = 0; j < M; ++j)
          std::cout << " " << buffer_iterators[0]->data()[i * buf_size + j];
        std::cout << std::endl;
        //      buffer_iterators[0]->data()[i*buf_size ]+=14;
      }
    }
    if (yy_constant_stride){
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, M, N, K, 1, buffer_iterators[0]->data(), buf_size,
                  alphas.data().data(), N, 1, yy_data + container_offset, yy_stride);
    }

    else{
      for (size_t i=0; 0<yy.size(); ++i){
        cblas_dgemv(CblasColMajor, CblasTrans, M, K, 1, buffer_iterators[0]->data(), buf_size,
                      alphas.data().data()+(alphas.cols()*i), 1, 1, yy[i].get().local_buffer()->data(), 1);
      }
    }
    // non-uniform stride: 

    for (auto &iter : buffer_iterators)
      ++iter;
  }
  prof->stop();
}
/**
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
    DistrArray::value_type buffer[16384];
    BufferManager x_buf = BufferManager(xx.at(ii).get(), buffer, 16384);
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
*/
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
void gemm_outer_default(Handler &handler, const Matrix<typename Handler::value_type> alphas, const CVecRef<AR> &xx,
                        const VecRef<AL> &yy) {
  for (size_t ii = 0; ii < alphas.rows(); ++ii) {
    for (size_t jj = 0; jj < alphas.cols(); ++jj) {
      handler.axpy(alphas(ii, jj), xx.at(ii).get(), yy[jj].get());
    }
  }
}

template <class AL, class AR = AL>
Matrix<typename array::mapped_or_value_type_t<AL>> gemm_inner_distr_distr(const CVecRef<AL> &xx,
                                                                          const CVecRef<AR> &yy) {
  // const size_t spacing = 1;
  using value_type = typename array::mapped_or_value_type_t<AL>;
  auto mat = Matrix<value_type>({xx.size(), yy.size()});
  if (xx.size() == 0 || yy.size() == 0)
    return mat;
  auto prof = molpro::Profiler::single()->push("gemm_inner_distr_distr (unbuffered)");
  prof += mat.cols() * mat.rows() * xx.at(0).get().local_buffer()->size() * 2;
  for (size_t j = 0; j < mat.cols(); ++j) {
    auto loc_y = yy.at(j).get().local_buffer();
    for (size_t i = 0; i < mat.rows(); ++i) {
      auto loc_x = xx.at(i).get().local_buffer();
      mat(i, j) = std::inner_product(begin(*loc_x), end(*loc_x), begin(*loc_y), (value_type)0);
      // mat(i,j) += cblas_ddot(end(*loc_x) - begin(*loc_x), begin(*loc_x), spacing, begin(*loc_y), spacing);
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
  if (xx.size() == 0 || yy.size() == 0)
    return mat;
  auto prof = molpro::Profiler::single()->push("gemm_inner_distr_distr (buffered)");
  prof += mat.cols() * mat.rows() * xx.at(0).get().local_buffer()->size() * 2;

  DistrArray::value_type *buffer = new DistrArray::value_type[8192 * 2 * mat.cols()];
  std::vector<BufferManager> buffers;
  for (size_t j = 0; j < mat.cols(); ++j) {
    buffers.emplace_back(BufferManager(yy.at(j).get(), buffer + (8192 * j), 8192, BufferManager::buffertype::Double));
  }

  for (size_t j = 0; j < mat.cols(); ++j) {
    size_t offset = 0;
    for (auto buffer = buffers[j].begin(); buffer != buffers[j].end(); offset += buffers[j].chunk_size, ++buffer) {
      for (size_t i = 0; i < mat.rows(); ++i) {
        size_t buflen = buffer->cend() - buffer->cbegin();
        mat(i, j) +=
            cblas_ddot(buflen, buffer->cbegin(), spacing, xx.at(i).get().local_buffer()->data() + offset, spacing);
        // mat(i, j) += std::inner_product(buffer->cbegin(), buffer->cend(), xx.at(i).get().local_buffer()->data() +
        // offset, (value_type)0);
      }
    }
  }
#ifdef HAVE_MPI_H
  MPI_Allreduce(MPI_IN_PLACE, const_cast<value_type *>(mat.data().data()), mat.size(), MPI_DOUBLE, MPI_SUM,
                xx.at(0).get().communicator());
#endif
  delete[] buffer;
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
        mat(i, j) += std::inner_product(buffer->cbegin(), buffer->cend(), xx.at(i).get().local_buffer()->data() +
offset, (value_type)0);
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
  if (xx.size() == 0 || yy.size() == 0)
    return mat;
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
  MPI_Allreduce(MPI_IN_PLACE, const_cast<value_type *>(mat.data().data()), mat.size(), MPI_DOUBLE, MPI_SUM,
                xx.at(0).get().communicator());
#endif
  return mat;
}

template <class Handler, class AL, class AR = AL>
Matrix<typename Handler::value_type> gemm_inner_default(Handler &handler, const CVecRef<AL> &xx,
                                                        const CVecRef<AR> &yy) {
  auto mat = Matrix<typename Handler::value_type>({xx.size(), yy.size()});
  if (xx.size() == 0 || yy.size() == 0)
    return mat;
  for (size_t ii = 0; ii < mat.rows(); ++ii) {
    for (size_t jj = 0; jj < mat.cols(); ++jj) {
      mat(ii, jj) = handler.dot(xx.at(ii).get(), yy.at(jj).get());
    }
  }
  return mat;
}

} // namespace molpro::linalg::array::util

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_GEMM_H
