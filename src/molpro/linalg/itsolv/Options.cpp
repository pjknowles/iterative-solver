#include "Options.h"

namespace molpro::linalg::itsolv {

void Options::copy(const Options& options) {
  convergence_threshold = options.convergence_threshold;
  n_roots = options.n_roots;
}

Options::Options(const std::map<std::string, std::string>& opt) {
  auto present = [&opt](const auto& key) { return opt.find(key) != opt.end(); };
  if (auto key = "n_roots"; present(key)) {
    n_roots = std::stoi(opt.at(key));
  };
  if (auto key = "convergence_threshold"; present(key)) {
    convergence_threshold = std::stod(opt.at(key));
  }
}

} // namespace molpro::linalg::itsolv