#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_DAVIDSONOPTIONS_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_DAVIDSONOPTIONS_H

#include <molpro/linalg/itsolv/options_map.h>

#include <cstddef>
#include <optional>

namespace molpro::linalg::itsolv {

struct DavidsonOptions {
  DavidsonOptions() = default;
  DavidsonOptions(const options_map& opt);

  //! Every reset_D iterations, discard the D space and move the current solutions into the Q space. 0 means never.
  //! Default: 0.
  std::optional<std::size_t> reset_D;
  //! Maximum size of the Q space immediately after a D-space reset. Default: unlimited, but never more than
  //! max_size_qspace.
  std::optional<std::size_t> reset_D_max_Q_size;
  //! Maximum number of vectors kept in the Q space; older vectors are folded into the D space. Default: unlimited.
  //! A limit bounds memory and disk use, but slows convergence of hard problems, and without reset_D small limits can
  //! prevent convergence altogether, so set reset_D (e.g. 10) together with it. For linear equations with an indefinite
  //! matrix, a limited Q space can make the projected problem nearly singular, so that the solution diverges (#628).
  std::optional<std::size_t> max_size_qspace;
  //! Minimum number of vectors kept in the Q space when it is truncated. Default: 0.
  std::optional<std::size_t> min_size_qspace;
  std::optional<double> contrib_thresh;
  std::optional<double> norm_thresh;
  std::optional<double> svd_thresh;
  std::optional<bool> hermiticity;
};

} // namespace molpro::linalg::itsolv

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_DAVIDSONOPTIONS_H
