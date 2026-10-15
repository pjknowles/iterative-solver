#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_GEMM_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_GEMM_H
#ifdef HAVE_MPI_H
#include <mpi.h>
#endif
#include "BufferManager.h"
#include "gather_all.h"
#include <functional>
#include <iterator>
#include <future>
#include <iostream>
#include <molpro/Options.h>
#include <molpro/Profiler.h>
#include <molpro/cblas.h>
#include <molpro/linalg/array/DistrArrayDisk.h>
#include <molpro/linalg/array/type_traits.h>
#include <molpro/linalg/itsolv/subspace/Matrix.h>
#include <molpro/linalg/itsolv/wrap.h>
#include <molpro/linalg/itsolv/wrap_util.h>
#include <molpro/linalg/options.h>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace molpro::linalg::array::util {

using molpro::linalg::itsolv::CVecRef;
using molpro::linalg::itsolv::VecRef;
using molpro::linalg::itsolv::subspace::Matrix;

enum gemm_type { inner, outer };

//! Whether arrays of type A are held on disk, and so must be read through buffers rather than local_buffer()
template <class A>
inline constexpr bool is_disk_array_v = std::is_base_of_v<DistrArrayDisk, std::decay_t<A>>;

//! The section [lo, hi) of a distributed array held by this process
template <class A>
auto local_range(const A& a) {
  return a.distribution().range(BufferManager<std::decay_t<A>>::rank_in(a));
}

// Building blocks for the buffered operations

//! Returns the constant distance between consecutive arrays at pointers p, if they are the columns of a matrix with n
//! rows, or 0 if they are not. A single array is given the stride n.
template <class T>
size_t constant_stride(const std::vector<T*>& p, size_t n) {
  if (p.size() < 2)
    return n;
  const auto stride = p[1] - p[0];
  for (size_t j = 1; j + 1 < p.size(); ++j)
    if (p[j + 1] - p[j] != stride)
      return 0;
  return stride >= std::ptrdiff_t(n) and stride > 0 ? size_t(stride) : 0;
}

//! Pointers to elements [g, ...) of arrays whose local sections, starting at global index lo, are at bases
template <class T>
std::vector<T*> offset_pointers(const std::vector<T*>& bases, size_t lo, size_t g) {
  auto p = bases;
  for (auto& el : p)
    el += g - lo;
  return p;
}

//! Pointers to the start of the local sections of distributed arrays held in memory
template <class T, class A>
std::vector<T*> local_buffer_pointers(const std::vector<std::reference_wrapper<A>>& aa) {
  std::vector<T*> p;
  for (const auto& a : aa)
    p.push_back(const_cast<T*>(a.get().local_buffer()->data()));
  return p;
}

//! Checks that the distributed arrays aa all have the local section range, so that they can be processed together
template <class A, class Range>
void check_local_ranges(const std::vector<std::reference_wrapper<A>>& aa, const Range& range, const std::string& name) {
  for (const auto& a : aa)
    if (local_range(a.get()) != range)
      throw std::out_of_range(name + ": the distributed arrays are distributed differently");
}

/*!
 * @brief Pointers to the elements of contiguous arrays rr, offset to the start of the local section of the
 * distributed array q
 * @throws std::out_of_range if any of rr differs from q in size
 */
template <class T, class R, class A>
std::vector<T*> contiguous_pointers(const std::vector<std::reference_wrapper<R>>& rr, const A& q,
                                    const std::string& name) {
  const auto lo = local_range(q).first;
  std::vector<T*> p;
  for (const auto& r : rr) {
    if (size_t(std::size(r.get())) != size_t(q.size()))
      throw std::out_of_range(name + ": contiguous and distributed arrays are of different size: " +
                              std::to_string(std::size(r.get())) + " " + std::to_string(q.size()));
    p.push_back(const_cast<T*>(std::data(r.get())) + lo);
  }
  return p;
}

/*!
 * @brief One chunk of n elements of a gemm between the arrays xx, which are the nx columns of x with leading
 * dimension ldx, and the arrays yy at pointers y.
 *
 * inner: alphadata[i + k * nx] += x_i . y_k
 * outer: y_k += sum_i alphadata[i * ny + k] * x_i
 *
 * If yy are evenly spaced this is a single dgemm, otherwise one dgemv for each of yy.
 */
template <class T>
void gemm_chunk(gemm_type type, size_t n, const T* x, size_t ldx, size_t nx, T* alphadata, const std::vector<T*>& y) {
  const size_t ny = y.size();
  if (n == 0 or nx == 0 or ny == 0)
    return;
  const auto dims = std::to_string(nx) + ", " + std::to_string(ny) + ", " + std::to_string(n);
  const auto y_stride = constant_stride(y, n);
  if (type == gemm_type::outer) {
    if (y_stride > 0) {
      auto prof = molpro::Profiler::single()->push("gemm_outer: cblas_dgemm dimensions " + dims);
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, ny, nx, 1, x, ldx, alphadata, ny, 1, y[0], y_stride);
    } else {
      auto prof = molpro::Profiler::single()->push("gemm_outer: cblas_dgemv dimensions " + dims);
      for (size_t k = 0; k < ny; ++k)
        cblas_dgemv(CblasColMajor, CblasNoTrans, n, nx, 1, x, ldx, alphadata + k, ny, 1, y[k], 1);
    }
  } else {
    if (y_stride > 0) {
      auto prof = molpro::Profiler::single()->push("gemm_inner: cblas_dgemm dimensions " + dims);
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, nx, ny, n, 1, x, ldx, y[0], y_stride, 1, alphadata, nx);
    } else {
      auto prof = molpro::Profiler::single()->push("gemm_inner: cblas_dgemv dimensions " + dims);
      for (size_t k = 0; k < ny; ++k)
        cblas_dgemv(CblasColMajor, CblasTrans, n, nx, 1, x, ldx, y[k], 1, 1, alphadata + k * nx, 1);
    }
  }
}

/*!
 * @brief gemm between distributed arrays xx, read through a BufferManager, and arrays yy, which can be held anywhere.
 *
 * Only the local section of xx is processed, chunk by chunk. For each chunk of global indices [g, g+n), y_pointers(g,
 * n) returns pointers to the same elements of each of yy, and y_done(g, n) is called once the chunk is complete, e.g.
 * to write back yy updated by gemm_type::outer.
 *
 * @param alphadata see gemm_chunk(), with nx = xx.size() and ny the number of yy
 */
template <class AD, class T = typename AD::value_type>
void gemm_buffered(T* alphadata, const CVecRef<AD>& xx,
                   const std::type_identity_t<std::function<std::vector<T*>(size_t g, size_t n)>>& y_pointers,
                   const std::function<void(size_t g, size_t n)>& y_done, gemm_type type) {
  auto prof = molpro::Profiler::single()->push("gemm_buffered");
  if (xx.empty())
    return;
  const auto [lo, hi] = local_range(xx.front().get());
  if (hi <= lo)
    return;
  auto options = molpro::linalg::options();
  auto number_of_buffers = options->parameter("GEMM_BUFFERS", 2);
  // BufferManager's buffer_size parameter is the size of ONE chunk (it already multiplies
  // by number_of_buffers internally to size the total pool) -- do not multiply by
  // number_of_buffers again here, or each chunk silently grows number_of_buffers-fold.
  const int buf_size = std::min(int(hi - lo), options->parameter("GEMM_PAGESIZE", 8192));
  molpro::Profiler::single()->start("gemm: buffer setup");
  BufferManager buffer(xx, buf_size, number_of_buffers);
  molpro::Profiler::single()->stop("gemm: buffer setup");
  for (auto buffer_iterator = buffer.begin(); buffer_iterator != buffer.end(); ++buffer_iterator) {
    const size_t g = lo + buffer.buffer_offset();
    const size_t n = buffer.buffer_size();
    gemm_chunk(type, n, buffer_iterator->data(), buffer.buffer_stride(), xx.size(), alphadata, y_pointers(g, n));
    if (y_done)
      y_done(g, n);
  }
}

/*!
 * @brief The overlap matrix between arrays yy, accessed as described in gemm_buffered(), and distributed arrays xx,
 * which are read through a BufferManager. Collective.
 * @return matrix of dimensions (ny, xx.size())
 */
template <class AD, class T = typename AD::value_type>
Matrix<T>
gemm_inner_buffered(size_t ny, const CVecRef<AD>& xx,
                    const std::type_identity_t<std::function<std::vector<T*>(size_t g, size_t n)>>& y_pointers,
                    const std::function<void(size_t g, size_t n)>& y_done = {}) {
  auto prof = molpro::Profiler::single()->push("gemm_inner (buffered)");
  auto alphas = Matrix<T>({ny, xx.size()});
  alphas.fill(0);
  if (ny == 0 or xx.empty())
    return alphas;
  const auto [lo, hi] = local_range(xx.front().get());
  prof += xx.size() * ny * (hi - lo) * 2;
  auto alphadata = const_cast<T*>(alphas.data().data());
  gemm_buffered(alphadata, xx, y_pointers, y_done, gemm_type::inner);
#ifdef HAVE_MPI_H
  // reduce over the arrays' own communicator, which need not be the global one
  MPI_Allreduce(MPI_IN_PLACE, alphadata, alphas.size(), MPI_DOUBLE, MPI_SUM, xx.front().get().communicator());
#endif
  return alphas;
}

//! Checks the dimensions of the coefficient matrix of gemm_outer
template <class T>
void check_gemm_outer_dimensions(const Matrix<T>& alphas, size_t nx, size_t ny, const std::string& name) {
  if (alphas.rows() != nx)
    throw std::out_of_range(name + ": dimensions of xx and alphas are different: " + std::to_string(alphas.rows()) +
                            " " + std::to_string(nx));
  if (alphas.cols() != ny)
    throw std::out_of_range(name + ": dimensions of yy and alphas are different: " + std::to_string(alphas.cols()) +
                            " " + std::to_string(ny));
}

//! The transpose of an overlap matrix, which for the hermitian inner product is its conjugate transpose
template <class T>
Matrix<T> conjugate_transpose(const Matrix<T>& m) {
  Matrix<T> result({m.cols(), m.rows()});
  molpro::linalg::itsolv::subspace::conjugate_transpose_copy(result, m);
  return result;
}

/*!
 * @brief yy += alphas^T xx, for arrays xx held in memory and distributed arrays yy held anywhere.
 *
 * Only the local section of yy is processed, in pages. Pages of yy without a local buffer are read with get() and
 * written back with put(), so that each element is read and written once. x_bases are pointers to the starts of the
 * local sections of xx; unless they are evenly spaced, each page of them is first packed into one matrix, so that the
 * update is always a single dgemm.
 */
template <class AD, class T = typename AD::value_type>
void gemm_outer_to_distr(const Matrix<T>& alphas, const std::vector<T*>& x_bases, const VecRef<AD>& yy) {
  if (x_bases.empty() or yy.empty())
    return;
  check_gemm_outer_dimensions(alphas, x_bases.size(), yy.size(), "gemm_outer_to_distr");
  auto prof = molpro::Profiler::single()->push("gemm_outer_to_distr");
  const auto [lo, hi] = local_range(yy.front().get());
  check_local_ranges(yy, std::make_pair(lo, hi), "gemm_outer_to_distr");
  prof += x_bases.size() * yy.size() * (hi - lo) * 2;
  const auto nx = x_bases.size();
  const auto ny = yy.size();
  const bool yy_on_disk = yy.front().get().disk_page_size() > 0;
  const size_t page = yy_on_disk ? yy.front().get().disk_page_size()
                                 : size_t(molpro::linalg::options()->parameter("GEMM_PAGESIZE", 8192));
  std::vector<T*> y_bases;
  if (not yy_on_disk)
    y_bases = local_buffer_pointers<T>(yy);
  std::vector<T> y_page, x_page;
  for (size_t g = lo; g < hi; g += page) {
    const size_t n = std::min(page, hi - g);
    std::vector<T*> y;
    if (yy_on_disk) {
      y_page.resize(n * ny);
      for (size_t k = 0; k < ny; ++k) {
        y.push_back(y_page.data() + k * n);
        yy[k].get().get(g, g + n, y.back());
      }
    } else
      y = offset_pointers(y_bases, lo, g);
    const auto x = offset_pointers(x_bases, lo, g);
    auto x_stride = constant_stride(x, n);
    const T* x_data = x[0];
    if (x_stride == 0) {
      x_page.resize(n * nx);
      for (size_t i = 0; i < nx; ++i)
        std::copy_n(x[i], n, x_page.data() + i * n);
      x_data = x_page.data();
      x_stride = n;
    }
    gemm_chunk(gemm_type::outer, n, x_data, x_stride, nx, const_cast<T*>(alphas.data().data()), y);
    if (yy_on_disk)
      for (size_t k = 0; k < ny; ++k)
        yy[k].get().put(g, g + n, y[k]);
  }
}

// Buffered: the disk arrays xx are read in chunks through a BufferManager

//! The pointers to, and on completion the write-back of, chunks of distributed arrays yy for gemm_buffered()
template <class AL, class T = typename array::mapped_or_value_type_t<AL>>
struct DistrChunks {
  DistrChunks(const VecRef<AL>& yy, gemm_type type) : yy(yy), type(type) {
    if (not yy.empty() and yy.front().get().disk_page_size() == 0) {
      lo = local_range(yy.front().get()).first;
      bases = local_buffer_pointers<T>(yy);
    }
  }
  std::vector<T*> pointers(size_t g, size_t n) {
    if (not bases.empty())
      return offset_pointers(bases, lo, g);
    // on disk: each chunk is read into memory, and written back after an outer product
    chunk.resize(yy.size() * n);
    std::vector<T*> p;
    for (size_t j = 0; j < yy.size(); ++j) {
      p.push_back(chunk.data() + j * n);
      yy[j].get().get(g, g + n, p.back());
    }
    return p;
  }
  void done(size_t g, size_t n) {
    if (bases.empty() and type == gemm_type::outer)
      for (size_t j = 0; j < yy.size(); ++j)
        yy[j].get().put(g, g + n, chunk.data() + j * n);
  }
  std::function<std::vector<T*>(size_t, size_t)> pointers_function() {
    return [this](size_t g, size_t n) { return pointers(g, n); };
  }
  std::function<void(size_t, size_t)> done_function() {
    return [this](size_t g, size_t n) { done(g, n); };
  }
  const VecRef<AL>& yy;
  const gemm_type type;
  size_t lo = 0;
  std::vector<T*> bases;
  std::vector<T> chunk;
};

template <class AL, class AD, std::enable_if_t<is_disk_array_v<AD>, int> = 0>
Matrix<typename array::mapped_or_value_type_t<AL>> gemm_inner_distr_distr(const CVecRef<AL>& yy,
                                                                          const CVecRef<AD>& xx) {
  auto non_const_yy = molpro::linalg::itsolv::const_cast_wrap(yy);
  DistrChunks<AL> chunks(non_const_yy, gemm_type::inner);
  return gemm_inner_buffered(yy.size(), xx, chunks.pointers_function());
}

template <class AD, class AL, std::enable_if_t<is_disk_array_v<AD> and not is_disk_array_v<AL>, int> = 0>
Matrix<typename array::mapped_or_value_type_t<AL>> gemm_inner_distr_distr(const CVecRef<AD>& xx,
                                                                          const CVecRef<AL>& yy) {
  // gemm_inner is the hermitian inner product, so the two orders are conjugate transposes
  return conjugate_transpose(gemm_inner_distr_distr(yy, xx));
}

template <class AL, class AD, std::enable_if_t<is_disk_array_v<AD>, int> = 0>
void gemm_outer_distr_distr(const Matrix<typename array::mapped_or_value_type_t<AL>> alphas, const CVecRef<AD>& xx,
                            const VecRef<AL>& yy) {
  if (yy.empty() or xx.empty())
    return;
  check_gemm_outer_dimensions(alphas, xx.size(), yy.size(), "gemm_outer_distr_distr");
  auto prof = molpro::Profiler::single()->push("gemm_outer_distr_distr (buffered)");
  prof += xx.size() * yy.size() * (local_range(xx[0].get()).second - local_range(xx[0].get()).first) * 2;
  DistrChunks<AL> chunks(yy, gemm_type::outer);
  gemm_buffered(const_cast<typename array::mapped_or_value_type_t<AL>*>(alphas.data().data()), xx,
                chunks.pointers_function(), chunks.done_function(), gemm_type::outer);
}

// Contiguous arrays rr, e.g. std::vector, holding all elements on every process, with distributed arrays qq.
// Each process works on the local section of qq; rr updated by gemm_outer are then replicated over all processes.

//! The overlap matrix (rr.size(), qq.size()) between contiguous arrays rr and distributed arrays qq. Collective.
template <class R, class AD>
Matrix<typename AD::value_type> gemm_inner_contiguous_distr(const CVecRef<R>& rr, const CVecRef<AD>& qq) {
  using T = typename AD::value_type;
  if (rr.empty() or qq.empty())
    return Matrix<T>({rr.size(), qq.size()});
  const auto& q = qq.front().get();
  check_local_ranges(qq, local_range(q), "gemm_inner_contiguous_distr");
  const auto lo = local_range(q).first;
  const auto bases = contiguous_pointers<T>(rr, q, "gemm_inner_contiguous_distr");
  return gemm_inner_buffered(rr.size(), qq, [&bases, lo](size_t g, size_t) { return offset_pointers(bases, lo, g); });
}

//! The overlap matrix (qq.size(), rr.size()) between distributed arrays qq and contiguous arrays rr. Collective.
template <class AD, class R>
Matrix<typename AD::value_type> gemm_inner_distr_contiguous(const CVecRef<AD>& qq, const CVecRef<R>& rr) {
  return conjugate_transpose(gemm_inner_contiguous_distr(rr, qq));
}

//! rr += alphas^T qq for contiguous arrays rr and distributed arrays qq. Collective.
template <class AD, class R>
void gemm_outer_distr_contiguous(const Matrix<typename AD::value_type>& alphas, const CVecRef<AD>& qq,
                                 const VecRef<R>& rr) {
  using T = typename AD::value_type;
  if (rr.empty() or qq.empty())
    return;
  check_gemm_outer_dimensions(alphas, qq.size(), rr.size(), "gemm_outer_distr_contiguous");
  auto prof = molpro::Profiler::single()->push("gemm_outer_distr_contiguous");
  const auto& q = qq.front().get();
  check_local_ranges(qq, local_range(q), "gemm_outer_distr_contiguous");
  const auto lo = local_range(q).first;
  prof += qq.size() * rr.size() * (local_range(q).second - lo) * 2;
  const auto bases = contiguous_pointers<T>(rr, q, "gemm_outer_distr_contiguous");
  gemm_buffered(
      const_cast<T*>(alphas.data().data()), qq,
      [&bases, lo](size_t g, size_t) { return offset_pointers(bases, lo, g); }, {}, gemm_type::outer);
  for (const auto& r : rr)
    gather_all(q.distribution(), q.communicator(), std::data(r.get()));
}

//! qq += alphas^T rr for contiguous arrays rr and distributed arrays qq. Local.
template <class R, class AD>
void gemm_outer_contiguous_distr(const Matrix<typename AD::value_type>& alphas, const CVecRef<R>& rr,
                                 const VecRef<AD>& qq) {
  using T = typename AD::value_type;
  if (rr.empty() or qq.empty())
    return;
  gemm_outer_to_distr(alphas, contiguous_pointers<T>(rr, qq.front().get(), "gemm_outer_contiguous_distr"), qq);
}

// Without buffers, for arrays held in memory

template <class AL, class AR = AL, std::enable_if_t<not is_disk_array_v<AL> and not is_disk_array_v<AR>, int> = 0>
Matrix<typename array::mapped_or_value_type_t<AL>> gemm_inner_distr_distr(const CVecRef<AL>& xx,
                                                                          const CVecRef<AR>& yy) {
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

template <class AL, class AR = AL, std::enable_if_t<not is_disk_array_v<AR>, int> = 0>
void gemm_outer_distr_distr(const Matrix<typename array::mapped_or_value_type_t<AL>> alphas, const CVecRef<AR>& xx,
                            const VecRef<AL>& yy) {
  if constexpr (is_disk_array_v<AL>) {
    // the arrays xx are in memory, and so can be addressed in place, while yy are on disk
    if (xx.empty() or yy.empty())
      return;
    check_local_ranges(xx, local_range(yy.front().get()), "gemm_outer_distr_distr");
    gemm_outer_to_distr(alphas, local_buffer_pointers<typename array::mapped_or_value_type_t<AL>>(xx), yy);
    return;
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
  // combine the sparse contributions to each yy, so that it is updated in one pass even when it is held on disk
  for (size_t ii = 0; ii < alphas.cols(); ++ii) {
    auto combined = std::decay_t<AR>{};
    for (size_t jj = 0; jj < alphas.rows(); ++jj)
      for (const auto& [i, v] : xx.at(jj).get())
        combined[i] += alphas(jj, ii) * v;
    yy[ii].get().axpy(1, combined);
  }
}

template <class AL, class AR = AL>
Matrix<typename array::mapped_or_value_type_t<AL>> gemm_inner_distr_sparse(const CVecRef<AL>& xx,
                                                                           const CVecRef<AR>& yy) {
  using value_type = typename array::mapped_or_value_type_t<AL>;
  auto mat = Matrix<value_type>({xx.size(), yy.size()});
  if (xx.size() == 0 || yy.size() == 0)
    return mat;
  const auto [lo, hi] = local_range(xx.at(0).get());
  for (size_t i = 0; i < mat.rows(); ++i) {
    const auto& x = xx.at(i).get();
    // an array in memory is read in place; one on disk has only the elements needed read
    const bool on_disk = x.disk_page_size() > 0;
    decltype(x.local_buffer()) loc_x;
    if (not on_disk)
      loc_x = x.local_buffer();
    for (size_t j = 0; j < mat.cols(); ++j) {
      mat(i, j) = 0;
      for (auto it = yy.at(j).get().lower_bound(lo); it != yy.at(j).get().end() and it->first < hi; ++it) {
        value_type xk;
        if (on_disk)
          x.get(it->first, it->first + 1, &xk);
        else
          xk = (*loc_x)[it->first - lo];
        mat(i, j) += xk * it->second;
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
