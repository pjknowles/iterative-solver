#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_ITERATIVESOLVER_H_OPTIONS_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_ITERATIVESOLVER_H_OPTIONS_H
#include <optional>

namespace molpro::linalg::itsolv {
//! Access point for different options in iterative solvers
struct Options {
  virtual ~Options() = default;
  std::optional<double> convergence_threshold; //!< convergence threshold
  std::optional<int> n_roots;                  //!< number of roots to solve for

  void copy(const Options& options) {
    convergence_threshold = options.convergence_threshold;
    n_roots = options.n_roots;
  }
};

struct ILinearEigensystemOptions : Options {};
struct ILinearEquationsOptions : Options {};
struct INoneLinearEquationsOptions : Options {};
struct IOptimizeOptions : Options {};

} // namespace molpro::linalg::itsolv
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_ITERATIVESOLVER_H_OPTIONS_H
