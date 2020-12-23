#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_NONLINEAREQUATIONSOPTIONSDIIS_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_NONLINEAREQUATIONSOPTIONSDIIS_H
#include <molpro/linalg/itsolv/Options.h>
#include <string>
namespace molpro::linalg::itsolv {
/*!
 * @brief Allows setting and getting of options for NonLinearEquationsDIIS instance via IterativeSolver base class
 */
struct NonLinearEquationsDIISOptions : public INonLinearEquationsOptions {
  NonLinearEquationsDIISOptions() = default;
  NonLinearEquationsDIISOptions(const options_map& opt);
  std::optional<int> max_size_qspace;
  std::optional<double> norm_thresh;
  std::optional<double> svd_thresh;
};

} // namespace molpro::linalg::itsolv

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_NONLINEAREQUATIONSOPTIONSDIIS_H
