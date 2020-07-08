#ifndef GCI_SRC_MOLPRO_GCI_ARRAY_UTIL_H
#define GCI_SRC_MOLPRO_GCI_ARRAY_UTIL_H
#include <memory>

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

template <typename T, class Compare> struct CompareAbs {
  constexpr bool operator()(const T &lhs, const T &rhs) const { return Compare()(std::abs(lhs), std::abs(rhs)); }
};

} // namespace util
} // namespace array
} // namespace gci
} // namespace molpro

#endif // GCI_SRC_MOLPRO_GCI_ARRAY_UTIL_H
