#include "OptimizeBFGSOptions.h"
#include "util.h"

namespace molpro::linalg::itsolv {

OptimizeBFGSOptions::OptimizeBFGSOptions(const options_map& opt) : OptimizeOptions(opt) {
  auto facet = util::StringFacet{};
  auto opt_upper = util::capitalize_keys(opt, facet);
  if (auto key = facet.toupper("max_size_qspace"); opt_upper.count(key)) {
    max_size_qspace = std::stoi(opt_upper.at(key));
  }
  if (auto key = facet.toupper("norm_thresh"); opt_upper.count(key)) {
    norm_thresh = std::stod(opt_upper.at(key));
  }
  if (auto key = facet.toupper("Wolfe_1"); opt_upper.count(key)) {
    Wolfe_1 = std::stod(opt_upper.at(key));
  }
  if (auto key = facet.toupper("Wolfe_2"); opt_upper.count(key)) {
    Wolfe_2 = std::stod(opt_upper.at(key));
  }
  if (auto key = facet.toupper("linesearch_grow_factor"); opt_upper.count(key)) {
    linesearch_grow_factor = std::stod(opt_upper.at(key));
  }
  if (auto key = facet.toupper("linesearch_tolerance"); opt_upper.count(key)) {
    linesearch_tolerance = std::stod(opt_upper.at(key));
  }
  if (auto key = facet.toupper("svd_thresh"); opt_upper.count(key)) {
    svd_thresh = std::stod(opt_upper.at(key));
  }
  if (auto key = facet.toupper("quasinewton_maximum_step"); opt_upper.count(key)) {
    quasinewton_maximum_step = std::stod(opt_upper.at(key));
  }
  if (auto key = facet.toupper("strong_Wolfe"); opt_upper.count(key)) {
    strong_Wolfe = std::stoi(opt_upper.at(key)) != 0;
  }
}

} // namespace molpro::linalg::itsolv