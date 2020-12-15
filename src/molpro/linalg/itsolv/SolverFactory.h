#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SOLVERFACTORY_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SOLVERFACTORY_H
#include <memory>
#include <molpro/linalg/itsolv/IterativeSolver.h>

namespace molpro::linalg::itsolv {
/*!
 * @brief Factory for creating instances of specific solver implementations from the corresponding Options object
 * @note The factory has to be instantiated for specific R,Q and P vector types
 *
 * Instantiating the factory
 * -------------------------
 * Assuming user defined containers are declared in ``UserDefinedVectors.h`` than the user can include the following
 * explicit instantiation for the factory in a `.cpp` file,
 *
 * @code{cpp}
 * // "solver_factory_instantiation.cpp"
 * #include <molpro/linalg/itsolv/SolverFactory-implementation.h>
 * #include "UserDefinedVectors.h"
 *
 * template class molpro::linalg::itsolv::SolverFactory<UserDefinedRvector, UserDefinedQvector, UserDefinedPvector>;
 * @endcode
 *
 * Using the factory
 * -----------------
 * @code{cpp}
 * #include <molpro/linalg/itsolv/SolverFactory.h>
 * #include "UserDefinedVectors.h"
 * #include "read_user_options.h"
 * using molpro::linalg::itsolv;
 * using R = UserDefinedRVector;
 * using Q = UserDefinedQVector;
 * using P = UserDefinedPVector;
 *
 * int main(int argn, char** argv){
 *      std::shared_ptr<Options> options = read_user_options(argn, argv);
 *
 *      std::shared_ptr<IterativeSolver> solver = SolverFactory<R,Q,P>::create(*options)
 *      solver.solve();
 *
 *      // Or, if a solver type is fixed
 *      std::shared_ptr<ILinearEigensystem> eigensystem_solver =
 *          SolverFactory<R,Q,P>::create<ILinearEigensystem>(*options) solver.solve();
 *      eigensystem_solver.solve();
 *      user_defined_print(eigensystem_solver.eigenvalues());
 * }
 * @code{cpp}
 *
 */
template <class Rvector, class Qvector, class Pvector>
class SolverFactory {
public:
  using R = Rvector;
  using Q = Qvector;
  using P = Pvector;

  static std::shared_ptr<IterativeSolver<R, Q, P>> create(const Options& options);

  template <template <typename, typename, typename> class SolverType>
  static std::shared_ptr<SolverType<R, Q, P>> create(const Options& options);
};
} // namespace molpro::linalg::itsolv
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SOLVERFACTORY_H
