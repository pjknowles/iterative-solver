#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_OPTIONS_MAP_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_OPTIONS_MAP_H
#include <map>
#include <string>

namespace molpro::linalg::itsolv {

using options_map = std::map<std::string, std::string>;

namespace util {

//! Returns options map with capitalized keys @note SFacet should be util::StringFacet, it is templated to reduce header
//! bloat
template <class SFacet>
std::map<std::string, std::string> capitalize_keys(const options_map& options, const SFacet& facet) {
  auto opt_upper = decltype(options){};
  for (const auto& [key, val] : options) {
    opt_upper[facet.toupper(key)] = val;
  }
  return opt_upper;
}

} // namespace util
} // namespace molpro::linalg::itsolv
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_OPTIONS_MAP_H
