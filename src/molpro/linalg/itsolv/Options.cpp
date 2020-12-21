#include "Options.h"

namespace molpro::linalg::itsolv {

void Options::copy(const Options& options) {
  convergence_threshold = options.convergence_threshold;
  n_roots = options.n_roots;
}

} // namespace molpro::linalg::itsolv