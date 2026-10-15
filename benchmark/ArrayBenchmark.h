#ifndef LINEARALGEBRA_BENCHMARK_ARRAYBENCHMARK_H_
#define LINEARALGEBRA_BENCHMARK_ARRAYBENCHMARK_H_
#include <molpro/mpi.h>
#ifdef HAVE_MPI_H
#include <molpro/linalg/array/DistrArrayMPI3.h>
#include <molpro/linalg/array/DistrArraySpan.h>
#include <molpro/linalg/array/util/Distribution.h>
#endif
#ifdef LINEARALGEBRA_ARRAY_GA
#include <ga-mpi.h>
#include <molpro/linalg/array/DistrArrayGA.h>
#endif
#include <molpro/Profiler.h>
#include <molpro/linalg/array/ArrayHandler.h>
#include <molpro/linalg/array/ArrayHandlerDDisk.h>
#include <molpro/linalg/array/ArrayHandlerDistr.h>
#include <molpro/linalg/array/ArrayHandlerIterable.h>
#include <molpro/linalg/array/DistrArrayFile.h>
#include <molpro/linalg/array/DistrArraySpan.h>
#ifdef LINEARALGEBRA_ARRAY_HDF5
#include <molpro/linalg/array/DistrArrayHDF5.h>
#include <molpro/linalg/array/PHDF5Handle.h>
#include <molpro/linalg/array/util.h>
#include <molpro/linalg/array/util/temp_phdf5_handle.h>
#endif
#include "molpro/linalg/array/ArrayHandlerDDiskDistr.h"
#include "molpro/linalg/array/ArrayHandlerDistrDDisk.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <ostream>
#include <vector>
namespace molpro {
namespace linalg {
//using Fast = molpro::linalg::array::DistrArrayMPI3;
 using Fast = molpro::linalg::array::DistrArraySpan;
// using Fast = std::vector<double>;


template <class T>
auto allocate(size_t n) {
  auto prof = molpro::Profiler::single()->push("allocate()");
  auto ptr = std::make_unique<T>(n);
  ptr->fill(1.0);
  return ptr;
}
template <class T>
auto allocatev(size_t n, size_t nvec) {
  auto prof = molpro::Profiler::single()->push("allocatev()");
  std::vector<T> result;
  result.reserve(nvec);
  for (size_t i = 0; i < nvec; ++i) {
    result.emplace_back(n);
    result.back().fill(1.0);
  }
  return result;
}
// DistrArraySpan is genuinely distributed (its local buffer must match the
// same per-rank range that a co-operating Slow array such as DistrArrayMPI3
// uses), so the caller-supplied span must be sized to this rank's local
// chunk, not the global dimension.
size_t fast_local_size(size_t n) {
  auto distribution = molpro::linalg::array::util::make_distribution_spread_remainder<typename Fast::index_type>(
      n, molpro::mpi::size_global());
  auto range = distribution.range(molpro::mpi::rank_global());
  return range.second - range.first;
}
auto allocateFast(std::vector<typename Fast::value_type>& buffer,size_t n) {
  auto result = std::make_unique<Fast>(n, array::Span<typename Fast::value_type>(buffer.data(),fast_local_size(n)), molpro::mpi::comm_global());
  result->fill(1.0);
  return result;
}
auto allocatevFast(std::vector<typename Fast::value_type>& buffer, size_t n, size_t nvec) {
std::vector<Fast> result;
result.reserve(nvec);
  const auto local_n = fast_local_size(n);
  for (size_t i = 0; i < nvec; ++i) {
    result.emplace_back(n, molpro::linalg::array::Span<typename Fast::value_type>(&buffer.at(i*local_n),local_n), molpro::mpi::comm_global());
    result.back().fill(1.0);
  }
  return result;
}
#ifdef LINEARALGEBRA_ARRAY_MPI3
template <>
auto allocate<array::DistrArrayMPI3>(size_t n) {
  auto result = std::make_unique<array::DistrArrayMPI3>(n, molpro::mpi::comm_global());
  result->fill(1.0);
  return result;
}
template <>
auto allocatev<array::DistrArrayMPI3>(size_t n, size_t nvec) {
  std::vector<array::DistrArrayMPI3> result;
  result.reserve(nvec);
  for (size_t i = 0; i < nvec; ++i) {
    result.emplace_back(n, molpro::mpi::comm_global());
    result.back().fill(1.0);
  }
  return result;
}
#endif
#ifdef LINEARALGEBRA_ARRAY_GA
template <>
auto allocate<array::DistrArrayGA>(size_t n) {
  auto result = std::make_unique<array::DistrArrayGA>(n, GA_MPI_Comm());
  result->fill(1.0);
  return result;
}
template <>
auto allocatev<array::DistrArrayGA>(size_t n, size_t nvec) {
  std::vector<array::DistrArrayGA> result;
  result.reserve(nvec);
  for (size_t i = 0; i < nvec; ++i) {
    result.emplace_back(n, GA_MPI_Comm());
    result.back().fill(1.0);
  }
  return result;
}
#endif
#ifdef LINEARALGEBRA_ARRAY_HDF5
template <>
auto allocate<array::DistrArrayHDF5>(size_t n) {
  auto handle = std::make_shared<array::util::PHDF5Handle>(
      array::util::temp_phdf5_handle("benchmark", molpro::mpi::comm_global()));
  handle->open_file(array::util::HDF5Handle::Access::read_write);
  handle->open_group("/");
  handle->close_file();
  auto result = std::make_unique<array::DistrArrayHDF5>(handle, n);
  result->fill(1.0);
  return result;
}
template <>
auto allocatev<array::DistrArrayHDF5>(size_t n, size_t nvec) {
  std::vector<array::DistrArrayHDF5> result;
  result.reserve(nvec);
  for (size_t i = 0; i < nvec; ++i) {
    auto handle = std::make_shared<array::util::PHDF5Handle>(
        array::util::temp_phdf5_handle("benchmark", molpro::mpi::comm_global()));
    handle->open_file(array::util::HDF5Handle::Access::read_write);
    handle->open_group("/");
    handle->close_file();
    result.emplace_back(handle, n);
    result.back().fill(1.0);
  }
  return result;
}
#endif

/*!
 * @brief
 * @tparam Slow slow class
 */
template <class Slow = Fast>
class ArrayBenchmark {
public:
  std::string m_title;

private:
  size_t m_size;
  double m_target_seconds;
  std::vector<typename Fast::value_type> m_fast_buffer;
  std::unique_ptr<Slow> m_bufferSlow;
  std::unique_ptr<Fast> m_bufferFast;
  std::vector<Slow> m_buffervSlow;
  std::vector<Fast> m_buffervFast;
  molpro::Profiler& m_profiler;
  std::unique_ptr<array::ArrayHandler<Fast, Slow>> m_fast_slow_handler;
  std::unique_ptr<array::ArrayHandler<Slow, Fast>> m_slow_fast_handler;
  bool m_profile_individual;
  int m_mpi_size;

public:
  size_t m_repeat;
  explicit ArrayBenchmark(const std::string& title, std::unique_ptr<array::ArrayHandler<Fast, Slow>> fast_slow_handler,
                          std::unique_ptr<array::ArrayHandler<Slow, Fast>> slow_fast_handler, size_t n = 10000000,
                          size_t n_Slow = 1, size_t n_Fast = 1, bool profile_individual = false,
                          double target_seconds = 1)
      : m_title(title), m_size(n), m_target_seconds(target_seconds), m_fast_buffer(fast_local_size(n)*n_Fast), m_bufferSlow(allocate<Slow>(n)),
        m_bufferFast(allocateFast(m_fast_buffer, n)),  m_buffervSlow(allocatev<Slow>(n, n_Slow)),
        m_buffervFast(allocatevFast(m_fast_buffer, n, n_Fast)), m_profiler(*molpro::Profiler::single()),
        m_fast_slow_handler(std::move(fast_slow_handler)), m_slow_fast_handler(std::move(slow_fast_handler)),
        m_profile_individual(profile_individual), m_mpi_size(molpro::mpi::size_global()),
        m_repeat(std::max(1, int(1e9 * m_target_seconds / m_size))) {
    m_profiler.reset(title);
    m_profiler.set_max_depth(10);
  }

  void dot() {
    auto prof = m_profiler.push("dot");
    for (size_t i = 0; i < m_repeat; i++)
      m_fast_slow_handler->dot(*m_bufferFast, *m_bufferSlow);
    prof += m_size * m_repeat / m_mpi_size;
  }

  void gemm_inner() {
    auto prof = m_profiler.push("gemm_inner");
    const size_t rep = std::max(size_t(1), m_repeat / m_buffervFast.size() / m_buffervSlow.size());
    for (size_t i = 0; i < rep; i++)
      m_fast_slow_handler->gemm_inner(itsolv::cwrap(m_buffervFast), itsolv::cwrap(m_buffervSlow));
    //    prof += rep * m_size * m_buffervSlow.size() * m_buffervFast.size() / m_mpi_size;
  }

  void gemm_inner_transpose() {
    auto prof = m_profiler.push("gemm_inner (transpose)");
    const size_t rep = std::max(size_t(1), m_repeat / m_buffervFast.size() / m_buffervSlow.size());
    for (size_t i = 0; i < rep; i++)
          m_slow_fast_handler->gemm_inner(itsolv::cwrap(m_buffervSlow), itsolv::cwrap(m_buffervFast));
    //    prof += rep * m_size * m_buffervSlow.size() * m_buffervFast.size() / m_mpi_size;
  }

  void axpy() {
    auto prof = m_profiler.push("axpy");
    for (size_t i = 0; i < m_repeat; i++)
      m_fast_slow_handler->axpy(1.0, *m_bufferSlow, *m_bufferFast);
    prof += m_size * m_repeat / m_mpi_size;
  }

  void gemm_outer() {
    auto prof = m_profiler.push("gemm_outer");
    const size_t rep = std::max(size_t(1), m_repeat / m_buffervFast.size() / m_buffervSlow.size());
    itsolv::subspace::Matrix<double> alpha({m_buffervSlow.size(), m_buffervFast.size()});
    alpha.fill(1.0);
    for (size_t i = 0; i < rep; i++)
      m_fast_slow_handler->gemm_outer(alpha, itsolv::cwrap(m_buffervSlow), itsolv::wrap(m_buffervFast));
    //    prof += rep * m_size * m_buffervSlow.size() * m_buffervFast.size() / m_mpi_size;
  }


  void copy_construct() {
            auto prof = m_profiler.push("copy construct");
            //TODO is there an atomic copy constructor?
//            auto buffer = allocate<Slow>(m_bufferSlow->size());
            for (size_t i = 0; i < std::min(m_repeat,decltype(m_repeat)(1000)); i++)
//              decltype(*m_bufferSlow) copy(*m_bufferFast);
              *m_bufferSlow = m_slow_fast_handler->copy(*m_bufferFast);
            prof += m_size * m_repeat / m_mpi_size;
  }

  void copyout() {
    auto prof = m_profiler.push("copy Slow -> Fast");
    for (size_t i = 0; i < std::min(m_repeat,decltype(m_repeat)(1000)); i++)
      *m_bufferFast = m_fast_slow_handler->copy(*m_bufferSlow);
    prof += m_size * m_repeat / m_mpi_size;
  }

  void fill() {
    using scalar_type = typename Slow::value_type;
    auto prof = m_profiler.push("fill");
    for (size_t i = 0; i < m_repeat; i++)
      m_fast_slow_handler->fill(static_cast<scalar_type>(1), *m_bufferFast);
    prof += m_size * m_repeat / m_mpi_size;
  }

  void scal() {
    //    using scalar_type = typename Slow::value_type;
    if (m_profile_individual)
      for (size_t i = 0; i < m_repeat; i++) {
        auto prof = m_profiler.push("scal");
        m_fast_slow_handler->scal(1, *m_bufferFast);
        prof += m_size / m_mpi_size;
      }
    else {
      auto prof = m_profiler.push("scal");
      for (size_t i = 0; i < m_repeat; i++)
        m_fast_slow_handler->scal(1, *m_bufferFast);
      prof += m_size * m_repeat / m_mpi_size;
    }
  }

  void all() {
    fill();
    scal();
    copy_construct();
    // a DistrArraySpan, which is a view of memory supplied by the caller, cannot be constructed from a disk array
    if constexpr (not array::is_disk_v<Slow>)
      copyout();
    dot();
    // axpy(); // TODO make axpy work with DistrArraySpan
    gemm_inner();
    gemm_inner_transpose();
    gemm_outer();
  }
  molpro::Profiler& profiler() { return m_profiler; }

  //  friend std::ostream& operator<<(std::ostream& os, ArrayBenchmark<Slow, Fast>& obj);
};
template <class Slow>
std::ostream& operator<<(std::ostream& os, ArrayBenchmark<Slow>& obj) {
  os << obj.profiler();
  return os;
}

template <class Slow = Fast>
ArrayBenchmark<Slow> ArrayBenchmarkIterable(std::string title, size_t n = 10000000, size_t n_Slow = 1,
                                            size_t n_Fast = 1, bool profile_individual = false,
                                            double target_seconds = 1) {
  return ArrayBenchmark<Slow>(title, std::make_unique<array::ArrayHandlerIterable<Fast, Slow>>(),
                              std::make_unique<array::ArrayHandlerIterable<Slow, Fast>>(), n, n_Slow, n_Fast,
                              profile_individual, target_seconds);
}

#ifdef LINEARALGEBRA_ARRAY_MPI3
template <class Slow = array::DistrArrayMPI3>
ArrayBenchmark<Slow> ArrayBenchmarkDistributed(std::string title, size_t n = 10000000, size_t n_Slow = 1,
                                               size_t n_Fast = 1, bool profile_individual = false,
                                               double target_seconds = 1) {
  return ArrayBenchmark<Slow>(title, std::make_unique<array::ArrayHandlerDistr<Fast, Slow>>(),
                              std::make_unique<array::ArrayHandlerDistr<Slow, Fast>>(), n, n_Slow, n_Fast,
                              profile_individual, target_seconds);
}
#endif

template <class Slow = array::DistrArrayFile>
ArrayBenchmark<Slow> ArrayBenchmarkDDisk(std::string title, size_t n = 10000000, size_t n_Slow = 1, size_t n_Fast = 1,
                                         bool profile_individual = false, double target_seconds = 1) {
  return ArrayBenchmark<Slow>(title, std::make_unique<array::ArrayHandlerDistrDDisk<Fast, Slow>>(),
                              std::make_unique<array::ArrayHandlerDDiskDistr<Slow, Fast>>(), n, n_Slow, n_Fast,
                              profile_individual, target_seconds);
}

/*!
 * @brief Compares the cost of visiting every element of a DistrArray through the generic,
 * STL-compatible iterator interface (begin()/end(), one remote-memory-access call per element)
 * against DistrArray's own native interface (local_buffer() and fill(), which touch only
 * process-local memory and use a single collective call to combine results across processes).
 *
 * Reading is compared via std::accumulate, writing via std::fill/std::transform equivalents.
 * The iterator has no notion of "this process's share of the array", so a write through it has
 * to be done by a single process to stay race-free, whereas the native fill() operates on every
 * process's local share in parallel.
 */
template <class Array>
void benchmarkIteratorVsNative(const std::string& title, size_t n, double target_seconds = 1) {
  auto profiler = molpro::Profiler::single();
  profiler->reset(title);
  profiler->set_max_depth(10);
  auto arr = allocate<Array>(n);
  const auto mpi_size = molpro::mpi::size_global();
  // Per-element remote-memory-access through the iterator is dominated by fixed per-call
  // latency rather than element count, so (unlike the other benchmarks in this file) repeats
  // are small, fixed constants rather than scaled from an assumed throughput: at even a few
  // hundred repeats of a per-element RMA loop, runtime is already dominated by that latency.
  const size_t iterator_repeat = 3;
  const size_t native_repeat = std::min(size_t(200), std::max(size_t(10), size_t(200 * target_seconds)));

  double sum_iterator = 0;
  {
    auto prof = profiler->push("read: std::accumulate(begin, end) [iterator]");
    for (size_t i = 0; i < iterator_repeat; ++i)
      sum_iterator = std::accumulate(cbegin(*arr), cend(*arr), 0.0);
    prof += n * iterator_repeat / mpi_size;
  }

  double sum_native = 0;
  {
    auto prof = profiler->push("read: local_buffer() + MPI_Allreduce [native]");
    for (size_t i = 0; i < native_repeat; ++i) {
      auto buffer = arr->local_buffer();
      sum_native = std::accumulate(buffer->begin(), buffer->end(), 0.0);
#ifdef HAVE_MPI_H
      MPI_Allreduce(MPI_IN_PLACE, &sum_native, 1, MPI_DOUBLE, MPI_SUM, arr->communicator());
#endif
    }
    prof += n * native_repeat / mpi_size;
  }

  if (molpro::mpi::rank_global() == 0 && std::abs(sum_iterator - sum_native) > 1.0e-6 * double(n))
    std::cerr << title << ": iterator sum " << sum_iterator << " disagrees with native sum " << sum_native
              << std::endl;

  {
    auto prof = profiler->push("write: std::fill(begin, end) [iterator]");
    for (size_t i = 0; i < iterator_repeat; ++i) {
      if (molpro::mpi::rank_global() == 0)
        std::fill(begin(*arr), end(*arr), 2.0);
      arr->sync();
    }
    prof += n * iterator_repeat / mpi_size;
  }

  {
    auto prof = profiler->push("write: fill() [native]");
    for (size_t i = 0; i < native_repeat; ++i)
      arr->fill(3.0);
    prof += n * native_repeat / mpi_size;
  }

  if (molpro::mpi::rank_global() == 0)
    std::cout << *profiler;
}

} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_BENCHMARK_ARRAYBENCHMARK_H_
