#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_LINEAREIGENSYSTEMOPTIONS_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_LINEAREIGENSYSTEMOPTIONS_H
#include <molpro/linalg/itsolv/Options.h>

namespace molpro::linalg::itsolv {
/*!
 * @brief Allows setting and getting of options for LinearEigensystem instance via IterativeSolver base class
 *
 * Example
 * @code{.cpp}
 *   #include <molpro/linalg/itsolv/CastOptions.h>
 *   #include <molpro/linalg/itsolv/LinearEigensystemOptions.h>
 *   #include "create_solver.h"
 *   using molpro::linalg::itsolv;
 *
 *   // create_solver is implemented by the user with definition in .cpp to avoid recompilation
 *   std::shared_ptr<IterativeSolver<R,Q,P>> solver = create_solver(SolverType::LinearEigensystem);
 *
 *   // Setting implementation dependent options for the solver through its base class
 *   auto options = std::make_shared<LinearEigensystemOptions>();
 *   options->max_size_q = 10;
 *   options->hermiticity = true;
 *   options->svd_thresh = 1.0e-12;
 *   solver->set_options(options);
 *
 *   // Retrieving options
 *   std::shared_ptr<Options> options_itsolv =  solver->get_options();
 *   std::cout << "convergence threshold= " << options_itsolv->coonvergence_threshold.get();
 *   // To access options specific to LinearEigensistem we need to down-cast
 *   std::shared_ptr<LinearEigensystemOptions> options_lin_eig = CastOptions::LinearEigensystem(solver->get_options());
 *   std::cout << "svd threshold= " << options_p->svd_threshold.get() << std::endl;
 * @endcode
 */
struct LinearEigensystemOptions : public ILinearEigensystemOptions {
  std::optional<int> reset_D;
  std::optional<int> reset_D_max_Q_size;
  std::optional<int> max_size_qspace;
  std::optional<double> norm_thresh;
  std::optional<double> svd_thresh;
  std::optional<bool> hermiticity;
};

} // namespace molpro::linalg::itsolv

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_LINEAREIGENSYSTEMOPTIONS_H
