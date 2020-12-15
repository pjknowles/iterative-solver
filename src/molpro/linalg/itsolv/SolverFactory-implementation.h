#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SOLVERFACTORY_IMPLEMENTATION_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SOLVERFACTORY_IMPLEMENTATION_H
#include <molpro/linalg/itsolv/LinearEigensystem.h>
#include <molpro/linalg/itsolv/LinearEquations.h>
#include <molpro/linalg/itsolv/SolverFactory.h>

namespace molpro::linalg::itsolv {

template <class R, class Q, class P>
std::shared_ptr<IterativeSolver<R, Q, P>>
SolverFactory<IterativeSolver<R, Q, P>>::create(const Options& options,
                                                const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers) {
  Options* options_ptr = &const_cast<Options&>(options);
  if (dynamic_cast<ILinearEquationsOptions*>(options_ptr)) {
    return SolverFactory<ILinearEquations<R, Q, P>>::create(options, handlers);
  } else if (dynamic_cast<ILinearEigensystemOptions*>(options_ptr)) {
    return SolverFactory<ILinearEigensystem<R, Q, P>>::create(options, handlers);
  }
}

template <class R, class Q, class P>
std::shared_ptr<ILinearEigensystem<R, Q, P>>
SolverFactory<ILinearEigensystem<R, Q, P>>::create(const Options& options,
                                                   const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers) {
  auto solver = std::make_shared<ILinearEigensystem<R, Q, P>>(handlers);
  solver->set_options(options);
  return solver;
}

template <class R, class Q, class P>
std::shared_ptr<ILinearEquations<R, Q, P>>
SolverFactory<ILinearEquations<R, Q, P>>::create(const Options& options,
                                                 const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers) {
  auto solver = std::make_shared<ILinearEquations<R, Q, P>>(handlers);
  solver->set_options(options);
  return solver;
}

} // namespace molpro::linalg::itsolv

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SOLVERFACTORY_IMPLEMENTATION_H
