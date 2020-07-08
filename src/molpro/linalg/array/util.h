#ifndef GCI_SRC_MOLPRO_GCI_ARRAY_UTIL_H
#define GCI_SRC_MOLPRO_GCI_ARRAY_UTIL_H
#include <memory>
#include <vector>

namespace molpro {
class Profiler;

namespace gci {
namespace array {
namespace util {

//! Starts profiler on construction and stops it on destruction, if passed pointer is not null
struct ScopeProfiler {
  std::weak_ptr<molpro::Profiler> m_prof;
  std::string m_name;
  ScopeProfiler(const std::shared_ptr<molpro::Profiler> &prof, std::string name);
  ~ScopeProfiler();
  ScopeProfiler() = delete;
  ScopeProfiler(const ScopeProfiler &) = delete;
  ScopeProfiler &operator=(const ScopeProfiler &) = delete;
};

struct Distribution {
  Distribution(int n_proc, size_t dimension);
  //! Maps fist and last index in the array to a pair of processes encapsulating the corresponding buffer range
  std::pair<int, int> process_map(unsigned long lo, unsigned long hi);
  //! start and size for section of array local to each process
  std::vector<std::pair<unsigned long, size_t>> proc_buffer;
  size_t dim;
};
} // namespace util
} // namespace array
} // namespace gci
} // namespace molpro

#endif // GCI_SRC_MOLPRO_GCI_ARRAY_UTIL_H
