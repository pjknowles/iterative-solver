#include "create_solver.h"

#include <stdexcept>

#include <molpro/linalg/itsolv/LinearEigensystem.h>
#include <molpro/linalg/itsolv/LinearEquations.h>
#include <molpro/linalg/itsolv/Optimize.h>
#include <molpro/linalg/itsolv/IterativeSolver.h>



#include <molpro/linalg/itsolv/SolverFactory-implementation.h>

namespace molpro::linalg::itsolv {
template class SolverFactory<Rvector, Qvector, Pvector>;
}

namespace molpro::test {

std::pair<std::shared_ptr<ILinearEigensystem<Rvector, Qvector, Pvector>>,
          std::shared_ptr<molpro::linalg::itsolv::Logger>>
create_LinearEigensystem() {
  auto handlers = std::make_shared<molpro::linalg::itsolv::ArrayHandlers<Rvector, Qvector, Pvector>>();
  auto options = linalg::itsolv::LinearEigensystemOptions();
  auto factory = linalg::itsolv::SolverFactory<Rvector, Qvector, Pvector>{};
  auto solver = factory.create(options, handlers);
  auto solver_eigen =
      std::dynamic_pointer_cast<molpro::linalg::itsolv::LinearEigensystem<Rvector, Qvector, Pvector>>(solver);
  auto logger = solver_eigen->logger;
  return {solver, logger};
}

std::pair<std::shared_ptr<ILinearEquations<Rvector, Qvector, Pvector>>, std::shared_ptr<Logger>>
create_LinearEquations() {
  auto handlers = std::make_shared<molpro::linalg::itsolv::ArrayHandlers<Rvector, Qvector, Pvector>>();
  auto options = linalg::itsolv::LinearEquationsOptions();
  auto factory = linalg::itsolv::SolverFactory<Rvector, Qvector, Pvector>{};
  auto solver = factory.create(options, handlers);
  auto solver_eq =
      std::dynamic_pointer_cast<molpro::linalg::itsolv::LinearEquations<Rvector, Qvector, Pvector>>(solver);
  auto logger = solver_eq->logger;
  return {solver, logger};
}
std::pair<std::shared_ptr<IOptimize<Rvector, Qvector, Pvector>>, std::shared_ptr<Logger>>
create_Optimize() {
  auto handlers = std::make_shared<molpro::linalg::itsolv::ArrayHandlers<Rvector, Qvector, Pvector>>();
  auto solver = std::make_shared<molpro::linalg::itsolv::Optimize<Rvector, Qvector, Pvector>>(handlers);
  auto logger = solver->logger;
  return {solver, logger};
}

} // namespace molpro::test
