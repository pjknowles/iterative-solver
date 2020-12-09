#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_LINEAREQUATIONSOPTIONS_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_LINEAREQUATIONSOPTIONS_H
#include <molpro/linalg/itsolv/Options.h>
namespace molpro::linalg::itsolv {
/*!
 * @brief Allows setting and getting of options for Optimize instance via IterativeSolver base class
 */
struct OptimizeOptions : public Options {
  std::optional<int> max_size_qspace;
  std::optional<double> norm_thresh;
  std::optional<double> svd_thresh;
};

} // namespace molpro::linalg::itsolv

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_LINEAREQUATIONSOPTIONS_H
