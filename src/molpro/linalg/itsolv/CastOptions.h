#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_CASTOPTIONS_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_CASTOPTIONS_H
#include <memory>
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/LinearEigensystemOptions.h>
#include <stdexcept>

namespace molpro::linalg::itsolv {
//! Safely down-cast Options to one of the implementations
struct CastOptions {
private:
  template <class OptionsType>
  static std::shared_ptr<OptionsType> cast(const std::shared_ptr<Options>& options, std::string options_type_name) {
    auto lin_eig_options = std::shared_ptr<OptionsType>{};
    if (options) {
      lin_eig_options = std::dynamic_pointer_cast<OptionsType>(options);
      if (!lin_eig_options)
        throw std::runtime_error("Failed to cast Options to " + options_type_name);
    }
    return lin_eig_options;
  }

public:
  static std::shared_ptr<LinearEigensystemOptions> LinearEigensystem(const std::shared_ptr<Options>& options) {
    return cast<LinearEigensystemOptions>(options, "LinearEigensystemOptions");
  }
};
} // namespace molpro::linalg::itsolv

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_CASTOPTIONS_H
