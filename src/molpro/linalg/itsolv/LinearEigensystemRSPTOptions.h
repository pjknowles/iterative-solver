#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_LINEAREIGENSYSTEMRSPTOPTIONS_H_
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_LINEAREIGENSYSTEMRSPTOPTIONS_H_

#include <molpro/linalg/itsolv/Options.h>
namespace molpro::linalg::itsolv {
/*!
 * @brief Allows setting and getting of options for LinearEquationsRSPT instance via IterativeSolver base class
 */
struct LinearEquationsRSPTOptions : public LinearEquationsOptions {
  LinearEquationsRSPTOptions() = default;
  LinearEquationsRSPTOptions(const options_map& opt);

  std::optional<double> norm_thresh;
  std::optional<double> svd_thresh;
};

} // namespace molpro::linalg::itsolv

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_LINEAREIGENSYSTEMRSPTOPTIONS_H_
