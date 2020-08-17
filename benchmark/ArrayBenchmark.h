#ifndef LINEARALGEBRA_BENCHMARK_ARRAYBENCHMARK_H_
#define LINEARALGEBRA_BENCHMARK_ARRAYBENCHMARK_H_
#include <molpro/Profiler.h>
#include <molpro/linalg/array/ArrayHandler.h>
#include <molpro/linalg/array/ArrayHandlerDistr.h>
#include <molpro/linalg/array/ArrayHandlerIterable.h>
#include <molpro/linalg/array/DistrArrayMPI3.h>
#include <ostream>
#include <vector>
namespace molpro {
namespace linalg {

template <class T>
auto resize(T& obj, size_t n) -> decltype(obj.resize, void()) {
  obj.resize(n);
}
template <class T>
auto resize(T& obj, size_t n) -> decltype(obj, void()) {}

template <class L = std::vector<double>, class R = L>
class ArrayBenchmark {
  size_t m_size;
  double m_target_seconds;
  L m_bufferL;
  R m_bufferR;
  molpro::Profiler m_profiler;
  std::unique_ptr<array::ArrayHandler<L, R>> m_handler;

public:
  explicit ArrayBenchmark(std::string title, std::unique_ptr<array::ArrayHandler<L, R>> handler, size_t n = 10000000,
                          double target_seconds = 1)
      : m_size(n), m_target_seconds(target_seconds), m_profiler(title), m_handler(std::move(handler)) {
    //    m_profiler.push("ArrayBenchmark constructor"); // TODO this doesn't work properly
    m_profiler.start("ArrayBenchmark constructor");
    resize(m_bufferL, n);
    resize(m_bufferR, n);
    m_profiler.stop("ArrayBenchmark constructor");
  }

  void dot() {
    size_t repeat = std::max(1, int(1e9 * m_target_seconds / m_size));
    auto prof = m_profiler.push("dot");
    for (auto i = 0; i < repeat; i++)
      m_handler->dot(m_bufferL, m_bufferR);
    prof += 2 * m_size * repeat;
  }

  void axpy() {
    size_t repeat = std::max(1, int(1e9 * m_target_seconds / m_size));
    auto prof = m_profiler.push("axpy");
    for (auto i = 0; i < repeat; i++)
      m_handler->axpy(1.0, m_bufferL, m_bufferR);
    prof += 2 * m_size * repeat;
  }

  void all() {
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
ArrayBenchmark<L, R> ArrayBenchmarkII(std::string title, size_t n = 10000000, double target_seconds = 1) {
  return ArrayBenchmark<L, R>(title, std::make_unique<array::ArrayHandlerIterable<L, R>>(), n, target_seconds);
}

template <class L = std::vector<double>, class R = array::DistrArrayMPI3>
ArrayBenchmark<L, R> ArrayBenchmarkID(std::string title, size_t n = 10000000, double target_seconds = 1) {
  return ArrayBenchmark<L, R>(title, std::make_unique<array::ArrayHandlerDistr<L, R>>(), n, target_seconds);
}

template <class L = array::DistrArrayMPI3, class R = L>
ArrayBenchmark<L, R> ArrayBenchmarkDD(std::string title, size_t n = 10000000, double target_seconds = 1) {
  return ArrayBenchmark<L, R>(title, std::make_unique<array::ArrayHandlerDistr<L, R>>(), n, target_seconds);
}

} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_BENCHMARK_ARRAYBENCHMARK_H_
