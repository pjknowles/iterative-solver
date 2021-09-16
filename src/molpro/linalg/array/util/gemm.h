#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_GEMM_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_GEMM_H
#ifdef HAVE_MPI_H
#include <mpi.h>
#endif
#include <future>
#include <iostream>
#include <molpro/Options.h>
#include <molpro/Profiler.h>
#include <molpro/cblas.h>
#include <molpro/linalg/array/DistrArrayFile.h>
#include <molpro/linalg/array/type_traits.h>
#include <molpro/linalg/itsolv/subspace/Matrix.h>
#include <molpro/linalg/itsolv/wrap.h>
#include <molpro/linalg/itsolv/wrap_util.h>
#include <molpro/linalg/options.h>
#include <numeric>
#include <vector>

using molpro::linalg::itsolv::CVecRef;
using molpro::linalg::itsolv::VecRef;
using molpro::linalg::itsolv::subspace::Matrix;

namespace molpro::linalg::array::util {

enum gemm_type { inner, outer };

// Buffered

template <class AL>
Matrix<typename array::mapped_or_value_type_t<AL>> gemm_inner_distr_distr(const CVecRef<AL>& yy,
                                                                          const CVecRef<DistrArrayFile>& xx) {
  auto prof = molpro::Profiler::single()->push("gemm_inner_distr_distr (buffered)");
  if (not yy.empty())
    prof += xx.size() * yy.size() * yy[0].get().local_buffer()->size() * 2;
  using value_type = typename array::mapped_or_value_type_t<AL>;
  auto alphas = Matrix<value_type>({yy.size(), xx.size()});
  alphas.fill(0);
  auto non_const_yy = molpro::linalg::itsolv::const_cast_wrap(yy);
  auto alphadata = const_cast<value_type*>(alphas.data().data());
  gemm_distr_distr(alphadata, xx, non_const_yy, gemm_type::inner);
#ifdef HAVE_MPI_H
  MPI_Allreduce(MPI_IN_PLACE, alphadata, alphas.size(), MPI_DOUBLE, MPI_SUM, molpro::mpi::comm_global());
#endif
  return alphas;
}

template <class AL>
void gemm_outer_distr_distr(const Matrix<typename array::mapped_or_value_type_t<AL>> alphas,
                            const CVecRef<DistrArrayFile>& xx, const VecRef<AL>& yy) {
  if (yy.empty() or xx.empty())
    return;
  auto prof = molpro::Profiler::single()->push("gemm_outer_distr_distr (buffered)");
  if (not yy.empty())
    prof += xx.size() * yy.size() * yy[0].get().local_buffer()->size() * 2;
  if (alphas.rows() != xx.size())
    throw std::out_of_range(std::string{"gemm_outer_distr_distr: dimensions of xx and alphas are different: "} +
                            std::to_string(alphas.rows()) + " " + std::to_string(xx.size()));
  if (alphas.cols() != yy.size())
    throw std::out_of_range(std::string{"gemm_outer_distr_distr: dimensions of yy and alphas are different: "} +
                            std::to_string(alphas.cols()) + " " + std::to_string(yy.size()));
  gemm_distr_distr(const_cast<typename array::mapped_or_value_type_t<AL>*>(alphas.data().data()), xx, yy,
                   gemm_type::outer);
}

template <class AL>
void gemm_distr_distr(array::mapped_or_value_type_t<AL>* alphadata, const CVecRef<DistrArrayFile>& xx,
                      const VecRef<AL>& yy, gemm_type gemm_type) {

  if (xx.size() == 0 || yy.size() == 0) {
    return;
  }

  bool yy_constant_stride = true;
  int previous_stride = 0;
  int yy_stride = yy.front().get().local_buffer()->size();
  for (size_t j = 0; j < std::max((size_t)1, yy.size()) - 1; ++j) {
    auto unique_ptr_j = yy.at(j).get().local_buffer()->data();
    auto unique_ptr_jp1 = yy.at(j + 1).get().local_buffer()->data();
    yy_stride = unique_ptr_jp1 - unique_ptr_j;
    //    std::cout << "j="<<j<<" y_stride="<<y_stride<<std::endl;
    if (j > 0)
      yy_constant_stride = yy_constant_stride && (yy_stride == previous_stride);
    previous_stride = yy_stride;
  }
  yy_constant_stride = yy_constant_stride && (yy_stride > 0);

  auto options = molpro::linalg::options();
  const BufferManager::buffertype number_of_buffers = (options->parameter("GEMM_BUFFERS", 2) > 1)
                                                          ? BufferManager::buffertype::Double
                                                          : BufferManager::buffertype::Single;
  const int buf_size = options->parameter("GEMM_PAGESIZE", 8192) * number_of_buffers;
  //  std::cout << "buf_size=" << buf_size << " number_of_buffers=" << number_of_buffers << std::endl;

  std::vector<DistrArray::value_type> buffers_memory(buf_size * xx.size());
  std::vector<BufferManager> buffers;
  std::vector<BufferManager::Iterator> buffer_iterators;
  buffers.reserve(xx.size());
  buffer_iterators.reserve(xx.size());
  for (size_t j = 0; j < xx.size(); ++j) {
    buffers.emplace_back(BufferManager(xx.at(j).get(),
                                       Span<DistrArray::value_type>(buffers_memory.data() + (buf_size * j), buf_size),
                                       number_of_buffers));
    buffer_iterators.emplace_back(buffers[j].begin());
  }
  int current_buf_size; // = buffer_iterators.front()->size();
  const auto y_size = int(yy[0].get().local_buffer()->size());
  for (int container_offset = 0; container_offset < y_size; container_offset += current_buf_size) {
    current_buf_size = buffer_iterators.front()->size();
    if (gemm_type == gemm_type::outer) {
      if (yy_constant_stride and not yy.empty()) {
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, current_buf_size, yy.size(), xx.size(), 1,
                    buffer_iterators[0]->data(), buf_size, alphadata, yy.size(), 1,
                    yy[0].get().local_buffer()->data() + container_offset, yy_stride);
      } else { // non-uniform stride:
        for (size_t i = 0; i < yy.size(); ++i) {
          cblas_dgemv(CblasColMajor, CblasNoTrans, current_buf_size, xx.size(), 1, buffer_iterators[0]->data(),
                      buf_size, alphadata + i, yy.size(), 1, yy[i].get().local_buffer()->data() + container_offset, 1);
        }
      }
    } else if (gemm_type == gemm_type::inner) {
      if (yy_constant_stride and not yy.empty()) {
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, xx.size(), yy.size(), current_buf_size, 1,
                    buffer_iterators[0]->data(), buf_size, yy[0].get().local_buffer()->data() + container_offset,
                    yy_stride, 1, alphadata, xx.size());
      } else { // non-uniform stride:
        for (size_t k = 0; k < yy.size(); ++k) {
          cblas_dgemv(CblasColMajor, CblasTrans, current_buf_size, xx.size(), 1, buffer_iterators[0]->data(), buf_size,
                      yy[k].get().local_buffer()->data() + container_offset, 1, 1, alphadata + k * xx.size(), 1);
        }
      }
    }
    for (auto& iter : buffer_iterators)
      ++iter;
  }
}

// Without buffers

template <class AL, class AR = AL>
Matrix<typename array::mapped_or_value_type_t<AL>> gemm_inner_distr_distr(const CVecRef<AL>& xx,
                                                                          const CVecRef<AR>& yy) {
  if (std::is_same<AL, DistrArrayFile>::value){
    throw std::runtime_error("gemm_inner_distr_distr (unbuffered) called with DistrArrayFile (should never happen!)");
  }                                                            
  // const size_t spacing = 1;
  using value_type = typename array::mapped_or_value_type_t<AL>;
  auto mat = Matrix<value_type>({xx.size(), yy.size()});
  if (xx.size() == 0 || yy.size() == 0)
    return mat;
  auto prof = molpro::Profiler::single()->push("gemm_inner_distr_distr (unbuffered)");
  if (not xx.empty())
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
  MPI_Allreduce(MPI_IN_PLACE, const_cast<value_type*>(mat.data().data()), mat.size(), MPI_DOUBLE, MPI_SUM,
                xx.at(0).get().communicator());
#endif
  return mat;
}

template <class AL, class AR = AL>
void gemm_outer_distr_distr(const Matrix<typename array::mapped_or_value_type_t<AL>> alphas, const CVecRef<AR>& xx,
                            const VecRef<AL>& yy) {
  if (std::is_same<AL, DistrArrayFile>::value){
    throw std::runtime_error("gemm_outer_distr_distr (unbuffered) called with DistrArrayFile (should never happen!)");
  }
  auto prof = molpro::Profiler::single()->push("gemm_outer_distr_distr (unbuffered)");
  if (not yy.empty())
    prof += alphas.rows() * alphas.cols() * yy[0].get().local_buffer()->size() * 2;
  for (size_t ii = 0; ii < alphas.rows(); ++ii) {
    auto loc_x = xx.at(ii).get().local_buffer();
    for (size_t jj = 0; jj < alphas.cols(); ++jj) {
      auto loc_y = yy[jj].get().local_buffer();
      for (size_t i = 0; i < loc_y->size(); ++i)
        (*loc_y)[i] += alphas(ii, jj) * (*loc_x)[i];
    }
  }
}

// Sparse

template <class AL, class AR = AL>
void gemm_outer_distr_sparse(const Matrix<typename array::mapped_or_value_type_t<AL>> alphas, const CVecRef<AR>& xx,
                             const VecRef<AL>& yy) {
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

template <class AL, class AR = AL>
Matrix<typename array::mapped_or_value_type_t<AL>> gemm_inner_distr_sparse(const CVecRef<AL>& xx,
                                                                           const CVecRef<AR>& yy) {
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
  MPI_Allreduce(MPI_IN_PLACE, const_cast<value_type*>(mat.data().data()), mat.size(), MPI_DOUBLE, MPI_SUM,
                xx.at(0).get().communicator());
#endif
  return mat;
}

// Handlers

template <class Handler, class AL, class AR = AL>
void gemm_outer_default(Handler& handler, const Matrix<typename Handler::value_type> alphas, const CVecRef<AR>& xx,
                        const VecRef<AL>& yy) {
  for (size_t ii = 0; ii < alphas.rows(); ++ii) {
    for (size_t jj = 0; jj < alphas.cols(); ++jj) {
      handler.axpy(alphas(ii, jj), xx.at(ii).get(), yy[jj].get());
    }
  }
}

template <class Handler, class AL, class AR = AL>
Matrix<typename Handler::value_type> gemm_inner_default(Handler& handler, const CVecRef<AL>& xx,
                                                        const CVecRef<AR>& yy) {
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
