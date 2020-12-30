#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_LINEAREQUATIONSDAVIDSONOPTIONS_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_LINEAREQUATIONSDAVIDSONOPTIONS_H
#include <molpro/linalg/itsolv/Options.h>
namespace molpro::linalg::itsolv {
/*!
 * @brief Allows setting and getting of options for LinearEquationsDavidson instance via IterativeSolver base class
 */
struct LinearEquationsDavidsonOptions : public LinearEquationsOptions {
  LinearEquationsDavidsonOptions() = default;
  LinearEquationsDavidsonOptions(const options_map& opt);

  std::optional<int> reset_D;
  std::optional<int> reset_D_max_Q_size;
  std::optional<int> max_size_qspace;
  std::optional<double> norm_thresh;
  std::optional<double> svd_thresh;
  std::optional<bool> hermiticity;
  std::optional<double> augmented_hessian;
};

} // namespace molpro::linalg::itsolv

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_LINEAREQUATIONSDAVIDSONOPTIONS_H
