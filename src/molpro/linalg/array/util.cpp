#include "util.h"
#include <molpro/Profiler.h>

namespace molpro {
namespace linalg {
namespace array {
namespace util {

ScopeProfiler::ScopeProfiler(const std::shared_ptr<molpro::Profiler>& prof, std::string name)
    : m_prof(prof), m_name(std::move(name)) {
  if (auto p = m_prof.lock())
    p->start(m_name);
}
ScopeProfiler::~ScopeProfiler() {
  if (auto p = m_prof.lock())
    p->stop(m_name);
}

} // namespace util
} // namespace array
} // namespace linalg
} // namespace molpro
