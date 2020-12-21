#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_ITERATIVESOLVER_H_OPTIONS_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_ITERATIVESOLVER_H_OPTIONS_H
#include <optional>

namespace molpro::linalg::itsolv {
//! Access point for different options in iterative solvers
struct Options {
  virtual ~Options() = default;
  std::optional<double> convergence_threshold; //!< convergence threshold
  std::optional<int> n_roots;                  //!< number of roots to solve for

  /*!
   * @brief copies options from source object
   * @note this is only necessary for internal use and should not be implemented in subclasses.
   */
  void copy(const Options& source);
};

struct ILinearEigensystemOptions : Options {};
struct ILinearEquationsOptions : Options {};
struct INonLinearEquationsOptions : Options {};
struct IOptimizeOptions : Options {};

} // namespace molpro::linalg::itsolv
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_ITERATIVESOLVER_H_OPTIONS_H
