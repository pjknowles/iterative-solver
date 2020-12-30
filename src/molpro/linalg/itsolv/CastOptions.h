#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_CASTOPTIONS_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_CASTOPTIONS_H
#include <memory>
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/LinearEigensystemDavidsonOptions.h>
#include <molpro/linalg/itsolv/LinearEquationsDavidsonOptions.h>
#include <molpro/linalg/itsolv/NonLinearEquationsDIISOptions.h>
#include <molpro/linalg/itsolv/OptimizeBFGSOptions.h>
#include <molpro/linalg/itsolv/OptimizeSDOptions.h>
#include <stdexcept>

namespace molpro::linalg::itsolv {
//! Safely down-cast Options to one of the implementations
struct CastOptions {
private:
  template <class OptionsType>
  static std::shared_ptr<OptionsType> cast(const std::shared_ptr<Options>& options,
                                           const std::string& options_type_name) {
    auto lin_eig_options = std::shared_ptr<OptionsType>{};
    if (options) {
      lin_eig_options = std::dynamic_pointer_cast<OptionsType>(options);
      if (!lin_eig_options)
        throw std::runtime_error("Failed to cast Options to " + options_type_name);
    }
    return lin_eig_options;
  }

  template <class OptionsType>
  static const OptionsType& cast(const Options& options, const std::string& options_type_name) {
    try {
      return dynamic_cast<const OptionsType&>(options);
    } catch (std::bad_cast& err) {
      throw std::runtime_error("Failed to cast Options to " + options_type_name);
    }
  }

public:
  static std::shared_ptr<LinearEigensystemDavidsonOptions> LinearEigensystem(const std::shared_ptr<Options>& options) {
    return cast<LinearEigensystemDavidsonOptions>(options, "LinearEigensystemOptions");
  }
  static const LinearEigensystemDavidsonOptions& LinearEigensystem(const Options& options) {
    return cast<LinearEigensystemDavidsonOptions>(options, "LinearEigensystemOptions");
  }
  static LinearEigensystemDavidsonOptions& LinearEigensystem(Options& options) {
    const auto& opt = const_cast<const Options&>(options);
    return const_cast<LinearEigensystemDavidsonOptions&>(cast<LinearEigensystemDavidsonOptions>(options, "LinearEigensystemOptions"));
  }

  static std::shared_ptr<LinearEquationsDavidsonOptions> LinearEquations(const std::shared_ptr<Options>& options) {
    return cast<LinearEquationsDavidsonOptions>(options, "LinearEquationsOptions");
  }
  static const LinearEquationsDavidsonOptions& LinearEquations(const Options& options) {
    return cast<LinearEquationsDavidsonOptions>(options, "LinearEquationsOptions");
  }
  static LinearEquationsDavidsonOptions& LinearEquations(Options& options) {
    const auto& opt = const_cast<const Options&>(options);
    return const_cast<LinearEquationsDavidsonOptions&>(cast<LinearEquationsDavidsonOptions>(options, "LinearEquationsOptions"));
  }

  static std::shared_ptr<NonLinearEquationsDIISOptions>
  NonLinearEquationsDIIS(const std::shared_ptr<Options>& options) {
    return cast<NonLinearEquationsDIISOptions>(options, "NonLinearEquationsDIISOptions");
  }
  static const NonLinearEquationsDIISOptions& NonLinearEquationsDIIS(const Options& options) {
    return cast<NonLinearEquationsDIISOptions>(options, "NonLinearEquationsDIISOptions");
  }
  static NonLinearEquationsDIISOptions& NonLinearEquationsDIIS(Options& options) {
    const auto& opt = const_cast<const Options&>(options);
    return const_cast<NonLinearEquationsDIISOptions&>(cast<NonLinearEquationsDIISOptions>(options, "NonLinearEquationsDIISOptions"));
  }


  static std::shared_ptr<OptimizeBFGSOptions> OptimizeBFGS(const std::shared_ptr<Options>& options) {
    return cast<OptimizeBFGSOptions>(options, "OptimizeBFGSOptions");
  }
  static const OptimizeBFGSOptions& OptimizeBFGS(const Options& options) {
    return cast<OptimizeBFGSOptions>(options, "OptimizeBFGSOptions");
  }
  static OptimizeBFGSOptions& OptimizeBFGS(Options& options) {
    const auto& opt = const_cast<const Options&>(options);
    return const_cast<OptimizeBFGSOptions&>(cast<OptimizeBFGSOptions>(options, "OptimizeBFGSOptions"));
  }

  static std::shared_ptr<OptimizeSDOptions> OptimizeSD(const std::shared_ptr<Options>& options) {
    return cast<OptimizeSDOptions>(options, "OptimizeSDOptions");
  }
  static const OptimizeSDOptions& OptimizeSD(const Options& options) {
    return cast<OptimizeSDOptions>(options, "OptimizeSDOptions");
  }
  static OptimizeSDOptions& OptimizeSD(Options& options) {
    const auto& opt = const_cast<const Options&>(options);
    return const_cast<OptimizeSDOptions&>(cast<OptimizeSDOptions>(options, "OptimizeSDOptions"));
  }
};
} // namespace molpro::linalg::itsolv

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_CASTOPTIONS_H
