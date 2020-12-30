#include "OptimizeSDOptions.h"
#include "util.h"

namespace molpro::linalg::itsolv {

OptimizeSDOptions::OptimizeSDOptions(const options_map& opt) : OptimizeOptions(opt) {
  auto facet = util::StringFacet{};
  auto opt_upper = util::capitalize_keys(opt, facet);
}

} // namespace molpro::linalg::itsolv