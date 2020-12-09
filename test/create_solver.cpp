#include "create_solver.h"

#include <stdexcept>

#include <molpro/linalg/itsolv/LinearEigensystem.h>
#include <molpro/linalg/itsolv/LinearEquations.h>
#include <molpro/linalg/itsolv/Optimize.h>
#include <molpro/linalg/itsolv/IterativeSolver.h>




namespace molpro::test {

std::pair<std::shared_ptr<ILinearEigensystem<Rvector, Qvector, Pvector>>,
          std::shared_ptr<molpro::linalg::itsolv::Logger>>
create_LinearEigensystem() {
  auto handlers = std::make_shared<molpro::linalg::itsolv::ArrayHandlers<Rvector, Qvector, Pvector>>();
  auto solver = std::make_shared<molpro::linalg::itsolv::LinearEigensystem<Rvector, Qvector, Pvector>>(handlers);
  auto logger = solver->logger;
  return {solver, logger};
}
std::pair<std::shared_ptr<ILinearEquations<Rvector, Qvector, Pvector>>, std::shared_ptr<Logger>>
create_LinearEquations() {
  auto handlers = std::make_shared<molpro::linalg::itsolv::ArrayHandlers<Rvector, Qvector, Pvector>>();
  auto solver = std::make_shared<molpro::linalg::itsolv::LinearEquations<Rvector, Qvector, Pvector>>(handlers);
  auto logger = solver->logger;
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
