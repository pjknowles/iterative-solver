#include "util.h"
#include <molpro/Profiler.h>

namespace molpro::gci::array::util {

ScopeProfiler::ScopeProfiler(const std::shared_ptr<molpro::Profiler>& prof, std::string name)
    : m_prof(prof), m_name(std::move(name)) {
  if (auto p = m_prof.lock())
    p->start(m_name);
}
ScopeProfiler::~ScopeProfiler() {
  if (auto p = m_prof.lock())
    p->stop(m_name);
}

Distribution::Distribution(int n_proc, size_t dimension) : dim(dimension) {
  proc_buffer.reserve(n_proc);
  // First get even distribution
  auto block = dim / n_proc;
  // Spread the remainder over the first processes
  auto n_extra = dim % n_proc;
  for (size_t i = 0; i < n_extra; ++i)
    proc_buffer.emplace_back(i * (block + 1), block + 1);
  for (size_t i = n_extra; i < n_proc; ++i)
    proc_buffer.emplace_back(i * block, block);
}

std::pair<int, int> Distribution::process_map(unsigned long lo, unsigned long hi) {
  auto block = dim / proc_buffer.size();
  auto n_extra = dim % proc_buffer.size();
  auto size_large_chunk = (block + 1) * n_extra;
  auto find_process = [size_large_chunk, block](auto ind) -> int {
    if (ind < size_large_chunk)
      return (ind / (block + 1));
    else
      return ((ind - size_large_chunk) / (block));
  };
  return {find_process(lo), find_process(hi)};
}

} // namespace molpro::gci::array::util
