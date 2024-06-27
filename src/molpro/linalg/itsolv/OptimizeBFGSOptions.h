#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_OPTIMIZEOPTIONSBFGS_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_OPTIMIZEOPTIONSBFGS_H
#include <molpro/linalg/itsolv/Options.h>
#include <string>
namespace molpro::linalg::itsolv {
/*!
 * @brief Allows setting and getting of options for OptimizeBFGS instance via IterativeSolver base class
 */
struct OptimizeBFGSOptions : public OptimizeOptions {
  OptimizeBFGSOptions() = default;
  OptimizeBFGSOptions(const options_map& opt);
  std::optional<int> max_size_qspace;
  std::optional<double> norm_thresh;
  std::optional<double> svd_thresh;
  std::optional<bool> strong_Wolfe; //!< Whether to use strong or weak Wolfe conditions
  std::optional<double> Wolfe_1;    //!< Acceptance parameter for function value
  std::optional<double> Wolfe_2;    //!< Acceptance parameter for function gradient
  std::optional<double>
      linesearch_tolerance; //!< If the predicted line search is within tolerance of 1, don't bother taking it
  std::optional<double> linesearch_grow_factor; //!< If the predicted line search step is extrapolation, limit the step
                                                //!< to this factor times the current step
  std::optional<std::string> interpolant; //!< Name of interpolant along search line
};

} // namespace molpro::linalg::itsolv

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_OPTIMIZEOPTIONSBFGS_H
