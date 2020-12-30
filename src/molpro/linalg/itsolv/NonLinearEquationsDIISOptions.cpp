#include "NonLinearEquationsDIISOptions.h"
#include "util.h"

namespace molpro::linalg::itsolv {

NonLinearEquationsDIISOptions::NonLinearEquationsDIISOptions(const options_map& opt) : NonLinearEquationsOptions(opt) {
  auto facet = util::StringFacet{};
  auto opt_upper = util::capitalize_keys(opt, facet);
  if (auto key = facet.toupper("max_size_qspace"); opt_upper.count(key)) {
    max_size_qspace = std::stoi(opt_upper.at(key));
  }
  if (auto key = facet.toupper("norm_thresh"); opt_upper.count(key)) {
    norm_thresh = std::stod(opt_upper.at(key));
  }
  if (auto key = facet.toupper("svd_thresh"); opt_upper.count(key)) {
    svd_thresh = std::stod(opt_upper.at(key));
  }
}

} // namespace molpro::linalg::itsolv