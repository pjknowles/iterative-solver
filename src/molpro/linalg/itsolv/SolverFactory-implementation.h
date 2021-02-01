#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SOLVERFACTORY_IMPLEMENTATION_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SOLVERFACTORY_IMPLEMENTATION_H
#include "OptimizeSD.h"
#include "OptimizeSDOptions.h"
#include <molpro/linalg/itsolv/LinearEigensystemDavidson.h>
#include <molpro/linalg/itsolv/LinearEigensystemRSPT.h>
#include <molpro/linalg/itsolv/LinearEquationsDavidson.h>
#include <molpro/linalg/itsolv/NonLinearEquationsDIIS.h>
#include <molpro/linalg/itsolv/OptimizeBFGS.h>
#include <molpro/linalg/itsolv/SolverFactory.h>
#include <molpro/linalg/itsolv/subspace/SubspaceSolverDIIS.h>
#include <molpro/linalg/itsolv/subspace/SubspaceSolverOptBFGS.h>
#include <molpro/linalg/itsolv/subspace/SubspaceSolverOptSD.h>

namespace molpro::linalg::itsolv {

template <class R, class Q, class P>
std::unique_ptr<IterativeSolver<R, Q, P>>
SolverFactory<R, Q, P>::create(const Options& options, const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers) {
  Options* options_ptr = &const_cast<Options&>(options);
  if (auto options_child = dynamic_cast<LinearEigensystemOptions*>(options_ptr); options_child) {
    return create(*options_child, handlers);
  } else if (auto options_child = dynamic_cast<LinearEquationsOptions*>(options_ptr); options_child) {
    return create(*options_child, handlers);
  } else {
    return nullptr;
  }
}

template <class R, class Q, class P>
std::unique_ptr<LinearEigensystem<R, Q, P>>
SolverFactory<R, Q, P>::create(const LinearEigensystemOptions& options,
                               const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers) {
  Options* options_ptr = &const_cast<LinearEigensystemOptions&>(options);
  if (auto options_child = dynamic_cast<LinearEigensystemDavidsonOptions*>(options_ptr); options_child) {
    auto solver = std::make_unique<LinearEigensystemDavidson<R, Q, P>>(handlers);
    solver->set_options(options);
    return solver;
  }
  if (auto options_child = dynamic_cast<LinearEigensystemRSPTOptions*>(options_ptr); options_child) {
    auto solver = std::make_unique<LinearEigensystemRSPT<R, Q, P>>(handlers);
    solver->set_options(options);
    return solver;
  }
}

template <class R, class Q, class P>
std::unique_ptr<LinearEquations<R, Q, P>>
SolverFactory<R, Q, P>::create(const LinearEquationsOptions& options,
                               const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers) {
  auto solver = std::make_unique<LinearEquationsDavidson<R, Q, P>>(handlers);
  solver->set_options(options);
  return solver;
}

template <class R, class Q, class P>
std::unique_ptr<NonLinearEquations<R, Q, P>>
SolverFactory<R, Q, P>::create(const NonLinearEquationsOptions& options,
                               const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers) {
  Options* options_ptr = &const_cast<NonLinearEquationsOptions&>(options);
  if (auto options_child = dynamic_cast<NonLinearEquationsDIISOptions*>(options_ptr); options_child) {
    auto solver = std::make_unique<NonLinearEquationsDIIS<R, Q, P>>(handlers);
    solver->set_options(options);
    return solver;
  }
  throw std::runtime_error("Unimplemented solver method");
}

template <class R, class Q, class P>
std::unique_ptr<Optimize<R, Q, P>>
SolverFactory<R, Q, P>::create(const OptimizeOptions& options,
                               const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers) {
  Options* options_ptr = &const_cast<OptimizeOptions&>(options);
  if (auto options_child = dynamic_cast<OptimizeBFGSOptions*>(options_ptr); options_child) {
    auto solver = std::make_unique<OptimizeBFGS<R, Q, P>>(handlers);
    solver->set_options(options);
    return solver;
  }
  if (auto options_child = dynamic_cast<OptimizeSDOptions*>(options_ptr); options_child) {
    auto solver = std::make_unique<OptimizeSD<R, Q, P>>(handlers);
    solver->set_options(options);
    return solver;
  }
  throw std::runtime_error("Unimplemented solver method");
}

template <class R, class Q, class P>
std::unique_ptr<IterativeSolver<R, Q, P>>
SolverFactory<R, Q, P>::create(const std::string& method, const options_map& options,
                               const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers) {
  if (method == "LinearEigensystem") {
    return create(LinearEigensystemDavidsonOptions{options}, handlers);
  } else if (method == "LinearEquations") {
    return create(LinearEquationsDavidsonOptions{options}, handlers);
  } else {
    throw std::runtime_error("Method = " + method + ", is not implemented");
  }
}

} // namespace molpro::linalg::itsolv

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SOLVERFACTORY_IMPLEMENTATION_H
