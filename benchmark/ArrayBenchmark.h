#ifndef LINEARALGEBRA_BENCHMARK_ARRAYBENCHMARK_H_
#define LINEARALGEBRA_BENCHMARK_ARRAYBENCHMARK_H_
#include <molpro/Profiler.h>
#include <molpro/linalg/array/ArrayHandler.h>
#include <vector>
#include <ostream>
namespace molpro {
namespace linalg {

template <class L = std::vector<double>, class R = L>
class ArrayBenchmark {
  size_t m_size;
  double m_target_seconds;
  L m_bufferL;
  R m_bufferR;
  molpro::Profiler m_profiler;
  std::unique_ptr<array::ArrayHandler<L, R>> m_handler;

public:
  explicit ArrayBenchmark(std::string title, size_t n = 10000000, double target_seconds = 1)
      : m_size(n), m_target_seconds(target_seconds), m_profiler(title) {
//    m_profiler.push("ArrayBenchmark constructor"); // TODO this doesn't work properly
    m_profiler.start("ArrayBenchmark constructor");
    m_bufferL.resize(n);
    m_bufferR.resize(n);
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
  molpro::Profiler& profiler() { return m_profiler;}

//  friend std::ostream& operator<<(std::ostream& os, ArrayBenchmark<L, R>& obj);
};
template <class L, class R>
std::ostream& operator<<(std::ostream& os, ArrayBenchmark<L, R>& obj) {
  os << obj.profiler();
  return os;
}
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_BENCHMARK_ARRAYBENCHMARK_H_
