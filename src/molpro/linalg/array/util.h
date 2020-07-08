#ifndef GCI_SRC_MOLPRO_GCI_ARRAY_UTIL_H
#define GCI_SRC_MOLPRO_GCI_ARRAY_UTIL_H
#include <memory>
#include <vector>

namespace molpro {
class Profiler;

namespace linalg {
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

} // namespace util
} // namespace array
} // namespace linalg
} // namespace molpro

#endif // GCI_SRC_MOLPRO_GCI_ARRAY_UTIL_H
