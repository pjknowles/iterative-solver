#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SOLVERFACTORY_IMPLEMENTATION_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SOLVERFACTORY_IMPLEMENTATION_H
#include "OptimizeSD.h"
#include "OptimizeSDOptions.h"
#include <molpro/linalg/itsolv/LinearEigensystem.h>
#include <molpro/linalg/itsolv/LinearEquations.h>
#include <molpro/linalg/itsolv/NonLinearEquations.h>
#include <molpro/linalg/itsolv/OptimizeBFGS.h>
#include <molpro/linalg/itsolv/SolverFactory.h>
#include <molpro/linalg/itsolv/subspace/SubspaceSolverDIIS.h>
#include <molpro/linalg/itsolv/subspace/SubspaceSolverOptBFGS.h>
#include <molpro/linalg/itsolv/subspace/SubspaceSolverOptSD.h>

namespace molpro::linalg::itsolv {

template <class R, class Q, class P>
std::shared_ptr<IterativeSolver<R, Q, P>>
SolverFactory<R, Q, P>::create(const Options& options, const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers) {
  Options* options_ptr = &const_cast<Options&>(options);
  if (auto options_child = dynamic_cast<ILinearEigensystemOptions*>(options_ptr); options_child) {
    return create(*options_child, handlers);
  } else if (auto options_child = dynamic_cast<ILinearEquationsOptions*>(options_ptr); options_child) {
    return create(*options_child, handlers);
  } else {
    return nullptr;
  }
}

template <class R, class Q, class P>
std::shared_ptr<ILinearEigensystem<R, Q, P>>
SolverFactory<R, Q, P>::create(const ILinearEigensystemOptions& options,
                               const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers) {
  auto solver = std::make_shared<LinearEigensystem<R, Q, P>>(handlers);
  solver->set_options(options);
  return solver;
}

template <class R, class Q, class P>
std::shared_ptr<ILinearEquations<R, Q, P>>
SolverFactory<R, Q, P>::create(const ILinearEquationsOptions& options,
                               const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers) {
  auto solver = std::make_shared<LinearEquations<R, Q, P>>(handlers);
  solver->set_options(options);
  return solver;
}

template <class R, class Q, class P>
std::shared_ptr<INonLinearEquations<R, Q, P>>
SolverFactory<R, Q, P>::create(const INonLinearEquationsOptions& options,
                               const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers) {
  Options* options_ptr = &const_cast<INonLinearEquationsOptions&>(options);
  if (auto options_child = dynamic_cast<NonLinearEquationsDIISOptions*>(options_ptr); options_child) {
    auto solver = std::make_shared<NonLinearEquations<subspace::SubspaceSolverDIIS, R, Q, P>>(handlers);
    solver->set_options(options);
    return solver;
  }
  throw std::runtime_error("Unimplemented solver method");
}

template <class R, class Q, class P>
std::shared_ptr<IOptimize<R, Q, P>>
SolverFactory<R, Q, P>::create(const IOptimizeOptions& options,
                               const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers) {
  Options* options_ptr = &const_cast<IOptimizeOptions&>(options);
  if (auto options_child = dynamic_cast<OptimizeBFGSOptions*>(options_ptr); options_child) {
    auto solver = std::make_shared<OptimizeBFGS<R, Q, P>>(handlers);
    solver->set_options(options);
    return solver;
  }
  if (auto options_child = dynamic_cast<OptimizeSDOptions*>(options_ptr); options_child) {
    auto solver = std::make_shared<OptimizeSD<R, Q, P>>(handlers);
    solver->set_options(options);
    return solver;
  }
  throw std::runtime_error("Unimplemented solver method");
}

template <class R, class Q, class P>
std::shared_ptr<IterativeSolver<R, Q, P>>
SolverFactory<R, Q, P>::create(const std::string& method, const options_map& options,
                               const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers) {
  if (method == "LinearEigensystem") {
    return create(LinearEigensystemOptions{options}, handlers);
  } else if (method == "LinearEquations") {
    return create(LinearEquationsOptions{options}, handlers);
  } else {
    throw std::runtime_error("Method = " + method + ", is not implemented");
  }
}

} // namespace molpro::linalg::itsolv

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SOLVERFACTORY_IMPLEMENTATION_H
