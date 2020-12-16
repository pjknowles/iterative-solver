#ifndef LINEARALGEBRA_TEST_CREATE_SOLVER_H
#define LINEARALGEBRA_TEST_CREATE_SOLVER_H
#include <memory>

#include <molpro/linalg/itsolv/CastOptions.h>
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/Logger.h>
#include <molpro/linalg/itsolv/SolverFactory.h>

using Rvector = std::vector<double>;
using Qvector = std::vector<double>;
using Pvector = std::map<size_t, double>;
using molpro::linalg::itsolv::CastOptions;
using molpro::linalg::itsolv::ILinearEigensystem;
using molpro::linalg::itsolv::ILinearEquations;
using molpro::linalg::itsolv::IOptimize;
using molpro::linalg::itsolv::IterativeSolver;
using molpro::linalg::itsolv::Logger;

namespace molpro::test {

std::pair<std::shared_ptr<ILinearEigensystem<Rvector, Qvector, Pvector>>, std::shared_ptr<Logger>>
create_LinearEigensystem();

std::pair<std::shared_ptr<ILinearEquations<Rvector, Qvector, Pvector>>, std::shared_ptr<Logger>>
create_LinearEquations();

std::pair<std::shared_ptr<IOptimize<Rvector, Qvector, Pvector>>, std::shared_ptr<Logger>>
create_Optimize();

} // namespace molpro::test
#endif // LINEARALGEBRA_TEST_CREATE_SOLVER_H
