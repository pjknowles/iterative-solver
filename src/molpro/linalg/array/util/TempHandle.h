#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_TEMPHANDLE_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_TEMPHANDLE_H
#include <map>
#include <memory>
#include <string>

namespace molpro {
namespace linalg {
namespace array {
namespace util {

//! Static registrar of temporary handles
template <class HandleType>
struct TempHandle {
  static std::map<std::string, std::weak_ptr<HandleType>> handles;
};

template <class HandleType>
std::map<std::string, std::weak_ptr<HandleType>> TempHandle<HandleType>::handles;

} // namespace util
} // namespace array
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_UTIL_TEMPHANDLE_H
