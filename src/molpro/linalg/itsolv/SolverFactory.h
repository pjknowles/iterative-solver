#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SOLVERFACTORY_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SOLVERFACTORY_H
#include <memory>

#include "LinearEigensystemDavidsonOptions.h"
#include "LinearEigensystemRSPTOptions.h"
#include "LinearEquationsDavidsonOptions.h"
#include "OptimizeBFGSOptions.h"
#include "OptimizeSDOptions.h"
#include "util.h"
#include <molpro/linalg/itsolv/ArrayHandlers.h>
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/NonLinearEquationsDIISOptions.h>
#include <molpro/linalg/itsolv/options_map.h>

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
 * namespace molpro::linalg::itsolv {
 * template class SolverFactory<UserDefinedRvector, UserDefinedQvector, UserDefinedPvector>;
 * }
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
 *      std::shared_ptr<Options> options = read_user_options_general(argn, argv);
 *      auto handlers = std::make_shared<ArrayHandlers<R, Q, P>>();
 *      auto factory = SolverFactory<R, Q, P>{};
 *
 *      std::unique_ptr<IterativeSolver> solver = factory.create(*options, handlers);
 *      solver->solve();
 *
 *      // Or, if a specific solver type is needed
 *      std::shared_ptr<ILinearEigensystemOptions> eigen_options = read_user_options_eigen(argn, argv);
 *      std::shared_ptr<ILinearEigensystem> eigen_solver = factory.create(*eigen_options, handlers);
 *      eigen_solver.solve();
 *      user_defined_print(eigen_solver.eigenvalues());
 *
 *      // Or with options provided as a key1=value1,key2=value2,.. string, and using a free-function factory call, and
 * using default handlers where these support the given R, Q, P: auto eigen_solver = create_LinearEigensystem<R, Q,
 * P>("max_size_qspace=6, nroots=3, convergence_threshold=1e-4"); eigen_solver.solve();
 *
 *      // The factory functions for non-linear equations and optimisation take
 *      // an additional argument specifying the method:
 *      auto minimiser = create_Optimize<R, Q, P>("BFGS","convergence_threshold=1e-5");
 * }
 * @endcode
 *
 */

template <class R, class Q = R, class P = std::map<size_t, typename R::value_type>>
class SolverFactory {
public:
  virtual ~SolverFactory() = default;
  virtual std::unique_ptr<IterativeSolver<R, Q, P>> create(const Options& options,
                                                           const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers);

  virtual std::unique_ptr<LinearEigensystem<R, Q, P>>
  create(const LinearEigensystemOptions& options,
         const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers =
             std::make_shared<molpro::linalg::itsolv::ArrayHandlers<R, Q, P>>());

  virtual std::unique_ptr<LinearEquations<R, Q, P>>
  create(const LinearEquationsOptions& options, const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers =
                                                    std::make_shared<molpro::linalg::itsolv::ArrayHandlers<R, Q, P>>());

  virtual std::unique_ptr<NonLinearEquations<R, Q, P>>
  create(const NonLinearEquationsOptions& options = NonLinearEquationsOptions{},
         const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers =
             std::make_shared<molpro::linalg::itsolv::ArrayHandlers<R, Q, P>>());

  virtual std::unique_ptr<Optimize<R, Q, P>>
  create(const OptimizeOptions& options = OptimizeOptions{},
         const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers =
             std::make_shared<molpro::linalg::itsolv::ArrayHandlers<R, Q, P>>());

  virtual std::unique_ptr<IterativeSolver<R, Q, P>>
  create(const std::string& method, const options_map& options,
         const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers =
             std::make_shared<molpro::linalg::itsolv::ArrayHandlers<R, Q, P>>());
};

// free-function factory invocation
template <class R, class Q = R, class P = std::map<size_t, typename R::value_type>>
std::unique_ptr<LinearEigensystem<R, Q, P>>
create_LinearEigensystem(const LinearEigensystemOptions& options,
                         const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers =
                             std::make_shared<molpro::linalg::itsolv::ArrayHandlers<R, Q, P>>()) {
  return SolverFactory<R, Q, P>{}.create(options, handlers);
}

template <class R, class Q = R, class P = std::map<size_t, typename R::value_type>>
std::unique_ptr<LinearEigensystem<R, Q, P>>
create_LinearEigensystem(const std::string& method = "Davidson", const std::string& options = "",
                         const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers =
                             std::make_shared<molpro::linalg::itsolv::ArrayHandlers<R, Q, P>>()) {
  auto optionsmap = util::StringFacet::parse_keyval_string(options);
  if (method == "Davidson" or method.empty())
    return SolverFactory<R, Q, P>{}.create(LinearEigensystemDavidsonOptions{optionsmap}, handlers);
  if (method == "RSPT")
    return SolverFactory<R, Q, P>{}.create(LinearEigensystemRSPTOptions{optionsmap}, handlers);
  throw std::runtime_error("Unimplemented method " + method);
}

template <class R, class Q = R, class P = std::map<size_t, typename R::value_type>>
std::unique_ptr<LinearEquations<R, Q, P>>
create_LinearEquations(const LinearEquationsOptions& options,
                       const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers =
                           std::make_shared<molpro::linalg::itsolv::ArrayHandlers<R, Q, P>>()) {
  return SolverFactory<R, Q, P>{}.create(options, handlers);
}

template <class R, class Q = R, class P = std::map<size_t, typename R::value_type>>
std::unique_ptr<LinearEquations<R, Q, P>>
create_LinearEquations(const std::string& method = "Davidson", const std::string& options = "",
                       const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers =
                           std::make_shared<molpro::linalg::itsolv::ArrayHandlers<R, Q, P>>()) {
  auto optionsmap = util::StringFacet::parse_keyval_string(options);
  if (method == "Davidson" or method.empty())
    return SolverFactory<R, Q, P>{}.create(LinearEquationsDavidsonOptions{optionsmap}, handlers);
  throw std::runtime_error("Unimplemented method " + method);
}

template <class R, class Q = R, class P = std::map<size_t, typename R::value_type>>
std::tuple<bool, std::unique_ptr<LinearEquations<R, Q, P>>>
Solve_LinearEquations(const VecRef<R>& parameters, const VecRef<R>& actions, const Problem<R>& problem,
                      int verbosity = 0, bool generate_initial_guess = true, const std::string& method = "Davidson",
                      const std::string& options = "",
                      const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers =
                          std::make_shared<molpro::linalg::itsolv::ArrayHandlers<R, Q, P>>()) {
  auto solver = create_LinearEquations<R, Q, P>(method, options, handlers);
  solver->set_verbosity(verbosity);
  auto success = solver->solve(parameters, actions, problem, generate_initial_guess);
//  if (success)
//    solver->solution(parameters, actions);
  return std::make_tuple<bool, std::unique_ptr<LinearEquations<R, Q, P>>>(std::move(success), std::move(solver));
}

template <class R, class Q = R, class P = std::map<size_t, typename R::value_type>>
std::tuple<bool, std::unique_ptr<LinearEquations<R, Q, P>>>
Solve_LinearEquations(std::vector<R>& parameters, std::vector<R>& actions, const Problem<R>& problem, int verbosity = 0,
                      bool generate_initial_guess = true, const std::string& method = "Davidson",
                      const std::string& options = "",
                      const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers =
                          std::make_shared<molpro::linalg::itsolv::ArrayHandlers<R, Q, P>>()) {
  return Solve_LinearEquations<R, Q, P>(wrap(parameters), wrap(actions), problem, verbosity, generate_initial_guess,
                                        method, options, handlers);
}

template <class R, class Q = R, class P = std::map<size_t, typename R::value_type>>
auto Solve_LinearEquations(R& parameters, R& actions, const Problem<R>& problem, int verbosity = 0,
                           bool generate_initial_guess = true, const std::string& method = "Davidson",
                           const std::string& options = "",
                           const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers =
                               std::make_shared<molpro::linalg::itsolv::ArrayHandlers<R, Q, P>>()) {
  return Solve_LinearEquations<R, Q, P>(std::vector<std::reference_wrapper<R>>{std::ref(parameters)},
                                        std::vector<std::reference_wrapper<R>>{std::ref(actions)}, problem, verbosity,
                                        generate_initial_guess, method, options, handlers);
}

template <class R, class Q = R, class P = std::map<size_t, typename R::value_type>>
std::unique_ptr<NonLinearEquations<R, Q, P>>
create_NonLinearEquations(const NonLinearEquationsOptions& options = NonLinearEquationsOptions{},
                          const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers =
                              std::make_shared<molpro::linalg::itsolv::ArrayHandlers<R, Q, P>>()) {
  return SolverFactory<R, Q, P>{}.create(options, handlers);
}

template <class R, class Q = R, class P = std::map<size_t, typename R::value_type>>
std::unique_ptr<NonLinearEquations<R, Q, P>>
create_NonLinearEquations(const std::string& method = "DIIS", const std::string& options = "",
                          const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers =
                              std::make_shared<molpro::linalg::itsolv::ArrayHandlers<R, Q, P>>()) {
  auto optionsmap = util::StringFacet::parse_keyval_string(options);
  if (method == "DIIS" or method.empty())
    return SolverFactory<R, Q, P>{}.create(NonLinearEquationsDIISOptions{optionsmap}, handlers);
  throw std::runtime_error("Unimplemented method " + method);
}

template <class R, class Q = R, class P = std::map<size_t, typename R::value_type>>
std::unique_ptr<Optimize<R, Q, P>>
create_Optimize(const OptimizeOptions& options = OptimizeOptions{},
                const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers =
                    std::make_shared<molpro::linalg::itsolv::ArrayHandlers<R, Q, P>>()) {
  return SolverFactory<R, Q, P>{}.create(options, handlers);
}

template <class R, class Q = R, class P = std::map<size_t, typename R::value_type>>
std::unique_ptr<Optimize<R, Q, P>>
create_Optimize(const std::string& method = "BFGS", const std::string& options = "",
                const std::shared_ptr<ArrayHandlers<R, Q, P>>& handlers =
                    std::make_shared<molpro::linalg::itsolv::ArrayHandlers<R, Q, P>>()) {
  auto optionsmap = util::StringFacet::parse_keyval_string(options);
  if (method == "BFGS" or method.empty())
    return SolverFactory<R, Q, P>{}.create(OptimizeBFGSOptions{optionsmap}, handlers);
  if (method == "SD")
    return SolverFactory<R, Q, P>{}.create(OptimizeSDOptions{optionsmap}, handlers);
  throw std::runtime_error("Unimplemented method " + method);
}

} // namespace molpro::linalg::itsolv
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SOLVERFACTORY_H
