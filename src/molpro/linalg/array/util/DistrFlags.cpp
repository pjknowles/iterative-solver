#include "DistrFlags.h"

#include <stdexcept>
namespace molpro {
namespace linalg {
namespace array {
namespace util {

int DistrFlags::Proxy::value() const {}

void DistrFlags::Proxy::validate() const {
  if (m_dflags_alive && *m_dflags_alive && !m_dflags.empty())
    return;
  else
    throw std::logic_error("overlying DistrFlags object is not valid");
}

} // namespace util
} // namespace array
} // namespace linalg
} // namespace molpro
