#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_ITERATIVESOLVER_H_OPTIONS_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_ITERATIVESOLVER_H_OPTIONS_H
#include <molpro/linalg/itsolv/options_map.h>

#include <map>
#include <optional>
#include <string>

namespace molpro::linalg::itsolv {
//! Access point for different options in iterative solvers
struct Options {
  virtual ~Options() = default;
  Options() = default;
  /*!
   * @brief Initialises options from key/value strings. Unrecognised options are ignored.
   * @param opt keys and values are  option name and corresponding value
   */
  Options(const options_map& opt);
  std::optional<double> convergence_threshold; //!< convergence threshold
  std::optional<int> n_roots;                  //!< number of roots to solve for

  /*!
   * @brief copies options from source object
   * @note this is only necessary for internal use and should not be implemented in subclasses.
   */
  void copy(const Options& source);
};

struct LinearEigensystemOptions : Options {
  LinearEigensystemOptions() = default;
  LinearEigensystemOptions(const options_map& opt) : Options(opt) {}
};

struct LinearEquationsOptions : Options {
  LinearEquationsOptions() = default;
  LinearEquationsOptions(const options_map& opt) : Options(opt) {}
};

struct NonLinearEquationsOptions : Options {
  NonLinearEquationsOptions() = default;
  NonLinearEquationsOptions(const options_map& opt) : Options(opt) {}
};

struct OptimizeOptions : Options {
  OptimizeOptions() = default;
  OptimizeOptions(const options_map& opt) : Options(opt) {}
};

} // namespace molpro::linalg::itsolv
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_ITERATIVESOLVER_H_OPTIONS_H
