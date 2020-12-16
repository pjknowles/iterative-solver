#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SOLVERFACTORY_IMPLEMENTATION_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SOLVERFACTORY_IMPLEMENTATION_H
#include <molpro/linalg/itsolv/LinearEigensystem.h>
#include <molpro/linalg/itsolv/LinearEquations.h>
#include <molpro/linalg/itsolv/SolverFactory.h>

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

} // namespace molpro::linalg::itsolv

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SOLVERFACTORY_IMPLEMENTATION_H
