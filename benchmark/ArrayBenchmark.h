#ifndef LINEARALGEBRA_BENCHMARK_ARRAYBENCHMARK_H_
#define LINEARALGEBRA_BENCHMARK_ARRAYBENCHMARK_H_
#include <molpro/mpi.h>
#ifdef HAVE_MPI_H
#include <molpro/linalg/array/DistrArrayMPI3.h>
#endif
#ifdef LINEARALGEBRA_ARRAY_GA
#include <ga-mpi.h>
#include <ga.h>
#include <macdecls.h>
#include <molpro/linalg/array/DistrArrayGA.h>
#endif
#include <molpro/Profiler.h>
#include <molpro/linalg/array/ArrayHandler.h>
#include <molpro/linalg/array/ArrayHandlerDDisk.h>
#include <molpro/linalg/array/ArrayHandlerDistr.h>
#include <molpro/linalg/array/ArrayHandlerIterable.h>
#include <molpro/linalg/array/DistrArraySpan.h>
#ifdef LINEARALGEBRA_ARRAY_HDF5
#include <molpro/linalg/array/DistrArrayHDF5.h>
#include <molpro/linalg/array/PHDF5Handle.h>
#include <molpro/linalg/array/util.h>
#include <molpro/linalg/array/util/temp_phdf5_handle.h>
#endif
#include "molpro/linalg/array/ArrayHandlerDDiskDistr.h"
#include "molpro/linalg/array/ArrayHandlerDistrDDisk.h"
#include <ostream>
#include <vector>
namespace molpro {
namespace linalg {
using Fast = molpro::linalg::array::DistrArrayMPI3;
// using Fast = molpro::linalg::array::DistrArraySpan;
// using Fast = std::vector<double>;

template <class T>
auto allocate(size_t n) {
  return std::make_unique<T>(n);
}
template <class T>
auto allocatev(size_t n, size_t nvec) {
  std::vector<T> result;
  result.reserve(nvec);
  for (size_t i = 0; i < nvec; ++i)
    result.emplace_back(n);
  return result;
}
#ifdef LINEARALGEBRA_ARRAY_MPI3
template <>
auto allocate<array::DistrArrayMPI3>(size_t n) {
  auto result = std::make_unique<array::DistrArrayMPI3>(n, molpro::mpi::comm_global());
  return result;
}
template <>
auto allocatev<array::DistrArrayMPI3>(size_t n, size_t nvec) {
  std::vector<array::DistrArrayMPI3> result;
  for (size_t i = 0; i < nvec; ++i)
    result.emplace_back(n, molpro::mpi::comm_global());
  return result;
}
#endif
#ifdef LINEARALGEBRA_ARRAY_GA
template <>
auto allocate<array::DistrArrayGA>(size_t n) {
  auto result = std::make_unique<array::DistrArrayGA>(n, GA_MPI_Comm());
  return result;
}
template <>
auto allocatev<array::DistrArrayGA>(size_t n, size_t nvec) {
  std::vector<array::DistrArrayGA> result;
  for (size_t i = 0; i < nvec; ++i)
    result.emplace_back(n, GA_MPI_Comm());
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
  return result;
}
template <>
auto allocatev<array::DistrArrayHDF5>(size_t n, size_t nvec) {
  std::vector<array::DistrArrayHDF5> result;
  for (size_t i = 0; i < nvec; ++i) {
    auto handle = std::make_shared<array::util::PHDF5Handle>(
        array::util::temp_phdf5_handle("benchmark", molpro::mpi::comm_global()));
    handle->open_file(array::util::HDF5Handle::Access::read_write);
    handle->open_group("/");
    handle->close_file();
    result.emplace_back(handle, n);
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
  std::unique_ptr<Slow> m_bufferSlow;
  std::unique_ptr<Fast> m_bufferFast;
  std::vector<Slow> m_buffervSlow;
  std::vector<Fast> m_buffervFast;
  molpro::Profiler& m_profiler;
  std::unique_ptr<array::ArrayHandler<Fast, Slow>> m_handler;
  bool m_profile_individual;
  int m_mpi_size;

public:
  size_t m_repeat;
  explicit ArrayBenchmark(const std::string& title, std::unique_ptr<array::ArrayHandler<Fast, Slow>> handler,
                          size_t n = 10000000, size_t n_Slow = 1, size_t n_Fast = 1, bool profile_individual = false,
                          double target_seconds = 1)
      : m_title(title), m_size(n), m_target_seconds(target_seconds), m_bufferSlow(allocate<Slow>(n)),
        m_bufferFast(allocate<Fast>(n)), m_buffervSlow(allocatev<Slow>(n, n_Slow)),
        m_buffervFast(allocatev<Fast>(n, n_Fast)), m_profiler(*molpro::Profiler::single()),
        m_handler(std::move(handler)), m_profile_individual(profile_individual), m_mpi_size(molpro::mpi::size_global()),
        m_repeat(std::max(1, int(1e9 * m_target_seconds / m_size))) {
    m_profiler.reset(title);
    m_profiler.set_max_depth(10);
  }

  void dot() {
    auto prof = m_profiler.push("dot");
    for (size_t i = 0; i < m_repeat; i++)
      m_handler->dot(*m_bufferFast, *m_bufferSlow);
    prof += m_size * m_repeat / m_mpi_size;
  }

  void gemm_inner() {
    auto prof = m_profiler.push("gemm_inner");
    const size_t rep = std::max(size_t(1), m_repeat / m_buffervFast.size() / m_buffervSlow.size());
    for (size_t i = 0; i < rep; i++)
      m_handler->gemm_inner(itsolv::cwrap(m_buffervFast), itsolv::cwrap(m_buffervSlow));
    prof += rep * m_size * m_buffervSlow.size() * m_buffervFast.size() / m_mpi_size;
  }

  void axpy() {
    auto prof = m_profiler.push("axpy");
    for (size_t i = 0; i < m_repeat; i++)
      m_handler->axpy(1.0, *m_bufferSlow, *m_bufferFast);
    prof += m_size * m_repeat / m_mpi_size;
  }

  void gemm_outer() {
    auto prof = m_profiler.push("gemm_outer");
    const size_t rep = std::max(size_t(1), m_repeat / m_buffervFast.size() / m_buffervSlow.size());
    Matrix<double> alpha({m_buffervSlow.size(), m_buffervFast.size()});
    alpha.fill(1.0);
    for (size_t i = 0; i < rep; i++)
      m_handler->gemm_outer(alpha, itsolv::cwrap(m_buffervSlow), itsolv::wrap(m_buffervFast));
    prof += rep * m_size * m_buffervSlow.size() * m_buffervFast.size() / m_mpi_size;
  }

  void copyin() {
//        auto prof = m_profiler.push("copy Fast -> Slow");
//        auto buffer = allocate<Slow>(m_bufferSlow->size());
//        for (size_t i = 0; i < m_repeat / m_mpi_size; i++)
//          *buffer = m_handler->copy(*m_bufferFast);
//        prof += m_size * m_repeat / m_mpi_size;
  }

  void copyout() {
    auto prof = m_profiler.push("copy Slow -> Fast");
    auto buffer = allocate<Fast>(m_bufferFast->size());
    for (size_t i = 0; i < m_repeat / m_mpi_size; i++)
      *buffer = m_handler->copy(*m_bufferSlow);
    prof += m_size * m_repeat / m_mpi_size;
  }

  void fill() {
    using scalar_type = typename Slow::value_type;
    auto prof = m_profiler.push("fill");
    for (size_t i = 0; i < m_repeat; i++)
      m_handler->fill(static_cast<scalar_type>(1), *m_bufferFast);
    prof += m_size * m_repeat / m_mpi_size;
  }

  void scal() {
    //    using scalar_type = typename Slow::value_type;
    if (m_profile_individual)
      for (size_t i = 0; i < m_repeat; i++) {
        auto prof = m_profiler.push("scal");
        m_handler->scal(1, *m_bufferFast);
        prof += m_size / m_mpi_size;
      }
    else {
      auto prof = m_profiler.push("scal");
      for (size_t i = 0; i < m_repeat; i++)
        m_handler->scal(1, *m_bufferFast);
      prof += m_size * m_repeat / m_mpi_size;
    }
  }

  void all() {
    fill();
    scal();
    copyin();
    copyout();
    dot();
    axpy();
    gemm_inner();
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
  return ArrayBenchmark<Slow>(title, std::make_unique<array::ArrayHandlerIterable<Fast, Slow>>(), n, n_Slow, n_Fast,
                              profile_individual, target_seconds);
}

#ifdef LINEARALGEBRA_ARRAY_MPI3
template <class Slow = array::DistrArrayMPI3>
ArrayBenchmark<Slow> ArrayBenchmarkDistributed(std::string title, size_t n = 10000000, size_t n_Slow = 1,
                                               size_t n_Fast = 1, bool profile_individual = false,
                                               double target_seconds = 1) {
  return ArrayBenchmark<Slow>(title, std::make_unique<array::ArrayHandlerDistr<Fast, Slow>>(), n, n_Slow, n_Fast,
                              profile_individual, target_seconds);
}
#endif

template <class Slow = array::DistrArrayFile>
ArrayBenchmark<Slow> ArrayBenchmarkDDisk(std::string title, size_t n = 10000000, size_t n_Slow = 1, size_t n_Fast = 1,
                                         bool profile_individual = false, double target_seconds = 1) {
  return ArrayBenchmark<Slow>(title, std::make_unique<array::ArrayHandlerDistrDDisk<Fast, Slow>>(), n, n_Slow, n_Fast,
                              profile_individual, target_seconds);
}

} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_BENCHMARK_ARRAYBENCHMARK_H_
