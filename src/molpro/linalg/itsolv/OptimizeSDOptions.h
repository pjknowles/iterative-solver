#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_OPTIMIZESDOPTIONS_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_OPTIMIZESDOPTIONS_H
#include <molpro/linalg/itsolv/Options.h>
#include <string>
namespace molpro::linalg::itsolv {
/*!
 * @brief Allows setting and getting of options for OptimizeBFGS instance via IterativeSolver base class
 */
struct OptimizeSDOptions : public OptimizeOptions {
  OptimizeSDOptions() = default;
  OptimizeSDOptions(const options_map& opt);
};

} // namespace molpro::linalg::itsolv

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_OPTIMIZESDOPTIONS_H
