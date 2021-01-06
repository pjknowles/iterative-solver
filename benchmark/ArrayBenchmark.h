#ifndef LINEARALGEBRA_BENCHMARK_ARRAYBENCHMARK_H_
#define LINEARALGEBRA_BENCHMARK_ARRAYBENCHMARK_H_
#ifdef HAVE_MPI_H
#include <molpro/linalg/array/DistrArrayMPI3.h>
#include <mpi.h>
#endif
#ifdef LINEARALGEBRA_ARRAY_GA
#undef MPI_COMM_WORLD
#define MPI_COMM_WORLD GA_MPI_Comm()
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
#ifdef LINEARALGEBRA_ARRAY_HDF5
#include <molpro/linalg/array/DistrArrayHDF5.h>
#include <molpro/linalg/array/PHDF5Handle.h>
#include <molpro/linalg/array/util.h>
#include <molpro/linalg/array/util/temp_phdf5_handle.h>
#endif
#include <ostream>
#include <vector>
namespace molpro {
namespace linalg {

template <class T>
auto allocate(size_t n) {
  return std::make_unique<T>(n);
}
#ifdef LINEARALGEBRA_ARRAY_MPI3
template <>
auto allocate<array::DistrArrayMPI3>(size_t n) {
  auto result = std::make_unique<array::DistrArrayMPI3>(n, MPI_COMM_WORLD);
  return result;
}
#endif
#ifdef LINEARALGEBRA_ARRAY_GA
template <>
auto allocate<array::DistrArrayGA>(size_t n) {
  auto result = std::make_unique<array::DistrArrayGA>(n, GA_MPI_Comm());
  return result;
}
#endif
#ifdef LINEARALGEBRA_ARRAY_HDF5
template <>
auto allocate<array::DistrArrayHDF5>(size_t n) {
  auto handle = std::make_shared<array::util::PHDF5Handle>(array::util::temp_phdf5_handle("benchmark", MPI_COMM_WORLD));
  handle->open_file(array::util::HDF5Handle::Access::read_write);
  handle->open_group("/");
  handle->close_file();
  auto result = std::make_unique<array::DistrArrayHDF5>(handle, n);
  result->open_access();
  return result;
}
#endif

template <class L = std::vector<double>, class R = L>
class ArrayBenchmark {
  size_t m_size;
  double m_target_seconds;
  std::unique_ptr<L> m_bufferL;
  std::unique_ptr<R> m_bufferR;
  molpro::Profiler m_profiler;
  std::unique_ptr<array::ArrayHandler<L, R>> m_handler;
  bool m_profile_individual;
  int m_mpi_size;

public:
  size_t m_repeat;
  explicit ArrayBenchmark(std::string title, std::unique_ptr<array::ArrayHandler<L, R>> handler, size_t n = 10000000,
                          bool profile_individual = false, double target_seconds = 1)
      : m_size(n), m_target_seconds(target_seconds), m_bufferL(allocate<L>(n)), m_bufferR(allocate<R>(n)),
        m_profiler(title), m_handler(std::move(handler)), m_profile_individual(profile_individual), m_mpi_size(1),
        m_repeat(std::max(1, int(1e9 * m_target_seconds / m_size))) {
#ifdef HAVE_MPI_H
    MPI_Comm_size(MPI_COMM_WORLD, &m_mpi_size);
#endif
  }

  void dot() {
    auto prof = m_profiler.push("dot");
    for (auto i = 0; i < m_repeat; i++)
      m_handler->dot(*m_bufferL, *m_bufferR);
    prof += m_size * m_repeat / m_mpi_size;
  }

  void axpy() {
    auto prof = m_profiler.push("axpy");
    for (auto i = 0; i < m_repeat; i++)
      m_handler->axpy(1.0, *m_bufferL, *m_bufferR);
    prof += m_size * m_repeat / m_mpi_size;
  }

  void copy() {
    auto prof = m_profiler.push("copy");
    auto buffer = allocate<L>(m_bufferL->size());
    for (auto i = 0; i < m_repeat / m_mpi_size; i++)
      *buffer = m_handler->copy(*m_bufferR);
    prof += m_size * m_repeat / m_mpi_size;
  }

  void fill() {
    using scalar_type = typename L::value_type;
    auto prof = m_profiler.push("fill");
    for (auto i = 0; i < m_repeat; i++)
      m_handler->fill(static_cast<scalar_type>(1), *m_bufferL);
    prof += m_size * m_repeat / m_mpi_size;
  }

  void scal() {
    using scalar_type = typename L::value_type;
    if (m_profile_individual)
      for (auto i = 0; i < m_repeat; i++) {
        auto prof = m_profiler.push("scal");
        m_handler->scal(1, *m_bufferL);
        prof += m_size / m_mpi_size;
      }
    else {
      auto prof = m_profiler.push("scal");
      for (auto i = 0; i < m_repeat; i++)
        m_handler->scal(1, *m_bufferL);
      prof += m_size * m_repeat / m_mpi_size;
    }
  }

  void all() {
    fill();
    scal();
    copy();
    dot();
    axpy();
  }
  molpro::Profiler& profiler() { return m_profiler; }

  //  friend std::ostream& operator<<(std::ostream& os, ArrayBenchmark<L, R>& obj);
};
template <class L, class R>
std::ostream& operator<<(std::ostream& os, ArrayBenchmark<L, R>& obj) {
  os << obj.profiler();
  return os;
}

template <class L = std::vector<double>, class R = L>
ArrayBenchmark<L, R> ArrayBenchmarkIterable(std::string title, size_t n = 10000000, bool profile_individual = false,
                                            double target_seconds = 1) {
  return ArrayBenchmark<L, R>(title, std::make_unique<array::ArrayHandlerIterable<L, R>>(), n, profile_individual,
                              target_seconds);
}

#ifdef LINEARALGEBRA_ARRAY_MPI3
template <class L = array::DistrArrayMPI3, class R = L>
ArrayBenchmark<L, R> ArrayBenchmarkDistributed(std::string title, size_t n = 10000000, bool profile_individual = false,
                                               double target_seconds = 1) {
  return ArrayBenchmark<L, R>(title, std::make_unique<array::ArrayHandlerDistr<L, R>>(), n, profile_individual,
                              target_seconds);
}
#endif

#ifdef LINEARALGEBRA_ARRAY_HDF5
template <class L = array::DistrArrayHDF5, class R = L>
ArrayBenchmark<L, R> ArrayBenchmarkDDisk(std::string title, size_t n = 10000000, bool profile_individual = false,
                                         double target_seconds = 1) {
  return ArrayBenchmark<L, R>(title, std::make_unique<array::ArrayHandlerDDisk<L, R>>(), n, profile_individual,
                              target_seconds);
}
#endif

} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_BENCHMARK_ARRAYBENCHMARK_H_
