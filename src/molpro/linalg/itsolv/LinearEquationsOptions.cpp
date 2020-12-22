#include "LinearEquationsOptions.h"
#include "util.h"

namespace molpro::linalg::itsolv {

LinearEquationsOptions::LinearEquationsOptions(const options_map& opt) : ILinearEquationsOptions(opt) {
  auto facet = util::StringFacet{};
  auto opt_upper = util::capitalize_keys(opt, facet);
  if (auto key = facet.toupper("reset_D"); opt_upper.count(key)) {
    reset_D = std::stoi(opt_upper.at(key));
  };
  if (auto key = facet.toupper("reset_D_max_Q_size"); opt_upper.count(key)) {
    reset_D_max_Q_size = std::stoi(opt_upper.at(key));
  }
  if (auto key = facet.toupper("max_size_qspace"); opt_upper.count(key)) {
    max_size_qspace = std::stoi(opt_upper.at(key));
  }
  if (auto key = facet.toupper("norm_thresh"); opt_upper.count(key)) {
    norm_thresh = std::stod(opt_upper.at(key));
  }
  if (auto key = facet.toupper("svd_thresh"); opt_upper.count(key)) {
    svd_thresh = std::stod(opt_upper.at(key));
  }
  if (auto key = facet.toupper("hermiticity"); opt_upper.count(key)) {
    hermiticity = facet.tobool(opt_upper.at(key));
  }
  if (auto key = facet.toupper("augmented_hessian"); opt_upper.count(key)) {
    augmented_hessian = std::stod(opt_upper.at(key));
  }
}

} // namespace molpro::linalg::itsolv