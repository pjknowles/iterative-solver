#include "IterativeSolverC.h"
#include "molpro/Profiler.h"
#include <algorithm>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <vector>
#ifdef LINEARALGEBRA_ARRAY_GA
#include "ga-mpi.h"
#include "ga.h"
#endif
#include <molpro/mpi.h>

// LINEARALGEBRA_ARRAY_HDF5 and LINEARALGEBRA_Q_HDF5 are deliberately different
// symbols, and both matter here. The first says the HDF5 containers are built; the
// second says to use them for the Q space of the C, Fortran and Python bindings in
// place of flat files. Nothing sets the second, by design: HDF5 measured slower than
// flat files for storage that is private to the solver and never seen by the client,
// so the Q space uses files, and using HDF5 instead would need a CMake option of its
// own. DistrArrayFile is therefore included unconditionally, because that is what the
// selection below resolves to; putting it in an #else left it undeclared in precisely
// those builds that do have HDF5.
#include <molpro/linalg/array/DistrArrayFile.h>
#ifdef LINEARALGEBRA_ARRAY_HDF5
#include <molpro/linalg/array/DistrArrayHDF5.h>
#endif
#include <molpro/linalg/itsolv/Logger.h>
#include <molpro/linalg/array/DistrArrayMPI3.h>
#include <molpro/linalg/array/DistrArraySpan.h>
#include <molpro/linalg/array/Span.h>
#include <molpro/linalg/array/util/Distribution.h>
#include <molpro/linalg/array/util/gather_all.h>
#include <molpro/linalg/itsolv/LinearEigensystemDavidson.h>
#include <molpro/linalg/itsolv/LinearEquationsDavidson.h>
#include <molpro/linalg/itsolv/SolverFactory.h>
#include <molpro/linalg/itsolv/wrap.h>
#include <molpro/linalg/itsolv/Logger.h>

using molpro::Profiler;
using molpro::linalg::array::Span;
using molpro::linalg::array::util::Distribution;
using molpro::linalg::array::util::gather_all;
using molpro::linalg::array::util::make_distribution_spread_remainder;
using molpro::linalg::itsolv::ArrayHandlers;
using molpro::linalg::itsolv::cwrap;
using molpro::linalg::itsolv::IterativeSolver;
using molpro::linalg::itsolv::LinearEigensystem;
using molpro::linalg::itsolv::LinearEigensystemDavidson;
using molpro::linalg::itsolv::LinearEquations;
using molpro::linalg::itsolv::LinearEquationsDavidson;
using molpro::linalg::itsolv::NonLinearEquations;
using molpro::linalg::itsolv::Optimize;
using molpro::linalg::itsolv::wrap;

// using Rvector = molpro::linalg::array::DistrArrayMPI3;
using Rvector = molpro::linalg::array::DistrArraySpan;
// Not LINEARALGEBRA_ARRAY_HDF5: see the note on the includes above. Flat files are
// the faster choice here, so this resolves to DistrArrayFile unless someone opts in.
#ifdef LINEARALGEBRA_Q_HDF5
using Qvector = molpro::linalg::array::DistrArrayHDF5;
#else
using Qvector = molpro::linalg::array::DistrArrayFile;
#endif
using Pvector = std::map<size_t, double>;
// instantiate the factory
#include <molpro/linalg/itsolv/SolverFactory-implementation.h>
template class molpro::linalg::itsolv::SolverFactory<Rvector, Qvector, Pvector>;

using vectorP = std::vector<double>;
using molpro::linalg::itsolv::CVecRef;
using molpro::linalg::itsolv::VecRef;

// Solver instances are identified by integer handles. All C API functions act on the current instance, which is the
// most recently created one unless another has been chosen with IterativeSolverSelect. IterativeSolverFinalize
// removes the current instance and makes the most recently created remaining one current, so code that creates and
// finalizes solvers in last-in-first-out order (the Fortran interface) never needs to select explicitly.

namespace {
struct Instance {
  Instance(std::unique_ptr<IterativeSolver<Rvector, Qvector, Pvector>> solver, std::shared_ptr<Profiler> prof,
           size_t dimension, MPI_Comm comm)
      : solver(std::move(solver)), prof(std::move(prof)), dimension(dimension), comm(comm){};
  std::unique_ptr<IterativeSolver<Rvector, Qvector, Pvector>> solver;
  std::shared_ptr<Profiler> prof;
  apply_on_p_t apply_on_p_fort;
  size_t dimension;
  MPI_Comm comm;
  Distribution<size_t> distribution; //!< distribution of the vectors over the processes of comm
  std::unique_ptr<Qvector> diagonals;
  bool has_values = false;
  bool has_eigenvalues = false;
  //! Optimize only: maximise the objective, by minimising its negative. Values and gradients are negated on their way
  //! into and out of the solver, so that the caller always sees those of its own objective.
  bool maximise = false;
  int verbosity = 0; //!< print level requested by the caller
};
class InstanceRegistry {
public:
  bool empty() const { return m_instances.empty(); }
  //! The current instance. @warning Undefined if empty()
  Instance& top() { return m_instances.at(m_current); }
  int64_t current_handle() const { return empty() ? -1 : m_current; }
  //! Add an instance, which becomes current
  void emplace(Instance&& instance) {
    m_current = m_next_handle++;
    m_instances.emplace(m_current, std::move(instance));
    m_creation_order.push_back(m_current);
  }
  //! Make an instance current. @return false if there is no instance with this handle
  bool select(int64_t handle) {
    if (m_instances.count(handle) == 0)
      return false;
    m_current = handle;
    return true;
  }
  //! Remove an instance; if it was current, the most recently created remaining instance becomes current
  void erase(int64_t handle) {
    if (m_instances.erase(handle) == 0)
      return;
    m_creation_order.erase(std::find(m_creation_order.begin(), m_creation_order.end(), handle));
    if (m_current == handle)
      m_current = m_creation_order.empty() ? -1 : m_creation_order.back();
  }
  //! Remove the current instance
  void pop() { erase(m_current); }
  void clear() {
    m_instances.clear();
    m_creation_order.clear();
    m_current = -1;
  }

private:
  std::map<int64_t, Instance> m_instances;
  std::vector<int64_t> m_creation_order;
  int64_t m_current = -1;
  int64_t m_next_handle = 0; //!< handles are never reused, so a stale handle cannot select a different solver
};
InstanceRegistry instances;
void require_instance();
} // namespace

std::pair<size_t, size_t> DistrArrayDefaultRange() {
  auto& instance = instances.top();
  int mpi_rank;
  MPI_Comm_rank(instance.comm, &mpi_rank);
  return instance.distribution.range(mpi_rank);
}

//! Marks a range that the caller did not specify
constexpr size_t unspecified_range = std::numeric_limits<size_t>::max();

/*!
 * @brief Set how the instance's vectors are distributed over the processes of its communicator.
 *
 * If every process supplies its local range (range_begin != unspecified_range), the ranges must cover
 * [0, dimension) contiguously in rank order, and they are used. If none does, the default distribution is used. On
 * return, range_begin and range_end hold this process's range. Collective over the instance's communicator.
 */
void set_distribution(Instance& instance, size_t* range_begin, size_t* range_end) {
  int mpi_rank, comm_size;
  MPI_Comm_rank(instance.comm, &mpi_rank);
  MPI_Comm_size(instance.comm, &comm_size);
  const unsigned long long mine[3] = {*range_begin != unspecified_range, *range_begin, *range_end};
  std::vector<unsigned long long> all(3 * comm_size);
  MPI_Allgather(mine, 3, MPI_UNSIGNED_LONG_LONG, all.data(), 3, MPI_UNSIGNED_LONG_LONG, instance.comm);
  int n_supplied = 0;
  for (int r = 0; r < comm_size; ++r)
    n_supplied += int(all[3 * r]);
  if (n_supplied == 0) {
    instance.distribution = make_distribution_spread_remainder<size_t>(instance.dimension, comm_size);
  } else {
    if (n_supplied != comm_size)
      throw std::invalid_argument("IterativeSolver: a range must be given on every process or on none");
    std::vector<size_t> borders{0};
    for (int r = 0; r < comm_size; ++r) {
      if (all[3 * r + 1] != borders.back() || all[3 * r + 2] < all[3 * r + 1])
        throw std::invalid_argument("IterativeSolver: the ranges of the processes must be contiguous, in rank order, "
                                    "starting at 0");
      borders.push_back(all[3 * r + 2]);
    }
    if (borders.back() != instance.dimension)
      throw std::invalid_argument("IterativeSolver: the ranges of the processes must cover the whole dimension");
    instance.distribution = Distribution<size_t>(borders);
  }
  std::tie(*range_begin, *range_end) = instance.distribution.range(mpi_rank);
}

std::vector<Rvector> CreateDistrArray(size_t nvec, double* data) {
  auto& instance = instances.top();
  MPI_Comm ccomm = instance.comm;
  int mpi_rank, comm_size;
  MPI_Comm_rank(instance.comm, &mpi_rank);
  MPI_Comm_size(instance.comm, &comm_size);
  const auto& distr = instance.distribution;
  auto range = distr.range(mpi_rank);
  auto rn = range.second - range.first;
  std::vector<Rvector> c;
  c.reserve(nvec);
  for (size_t ivec = 0; ivec < nvec; ivec++) {
    c.emplace_back(std::make_unique<Distribution<Rvector::index_type>>(distr),
                   Span<typename Rvector::value_type>(&data[ivec * instance.dimension + range.first], rn), ccomm);
  }
  return c;
}

std::vector<Rvector> CreateDistrArray(size_t nvec, const double* data) {
  auto& instance = instances.top();
  MPI_Comm ccomm = instance.comm;
  int mpi_rank, comm_size;
  MPI_Comm_rank(instance.comm, &mpi_rank);
  MPI_Comm_size(instance.comm, &comm_size);
  const auto& distr = instance.distribution;
  auto range = distr.range(mpi_rank);
  auto rn = range.second - range.first;
  std::vector<Rvector> c;
  c.reserve(nvec);
  for (size_t ivec = 0; ivec < nvec; ivec++) {
    c.emplace_back(
        std::make_unique<Distribution<Rvector::index_type>>(distr),
        Span<typename Rvector::value_type>(&const_cast<double*>(data)[ivec * instance.dimension + range.first], rn),
        ccomm);
  }
  return c;
}

void DistrArraySynchronize(size_t nvec, std::vector<Rvector>& c, double* data) {
  auto& instance = instances.top();
  MPI_Comm ccomm = instance.comm;
  for (size_t ivec = 0; ivec < nvec; ivec++) {
    gather_all(c[ivec].distribution(), ccomm, &data[ivec * instance.dimension]);
  }
}

std::pair<size_t, size_t> DistrArrayGetRange(Rvector& rvec) {
  auto& instance = instances.top();
  int mpi_rank;
  MPI_Comm_rank(instance.comm, &mpi_rank);
  auto range = rvec.distribution().range(mpi_rank);
  return range;
}

extern "C" void IterativeSolverRange(size_t* range_begin, size_t* range_end) {
  std::tie(*range_begin, *range_end) = DistrArrayDefaultRange();
}

void apply_on_p_c(const std::vector<vectorP>& pvectors, const CVecRef<Pvector>& pspace, const VecRef<Rvector>& action) {
  auto& instance = instances.top();
  std::vector<size_t> ranges;
  size_t update_size = pvectors.size();
  ranges.reserve(update_size * 2);
  for (size_t k = 0; k < update_size; ++k) {
    ranges.push_back(DistrArrayGetRange(action[k].get()).first);
    ranges.push_back(DistrArrayGetRange(action[k].get()).second);
  }
  std::vector<double> pvecs_to_send;
  for (size_t i = 0; i < update_size; i++) {
    for (auto j : pvectors[i]) {
      pvecs_to_send.push_back(j);
    }
  }
  // The callback indexes the actions as full-length arrays and fills only this process's range, so it needs the
  // address of the first element of the full array, not of this process's local section
  auto local = action.front().get().local_buffer();
  double* full_action = local->data() - ranges.front();
  instance.apply_on_p_fort(pvecs_to_send.data(), full_action, update_size, ranges.data());
}

extern "C" void IterativeSolverLinearEigensystemInitialize(size_t nQ, size_t nroot, size_t* range_begin,
                                                           size_t* range_end, double thresh, double thresh_value,
                                                           int hermitian, int verbosity, const char* fname,
                                                           int64_t fcomm, const char* algorithm, const char* options) {
  std::shared_ptr<Profiler> profiler = nullptr;
  profiler = molpro::Profiler::single("MainProfiler");
  std::string pname(fname);
  MPI_Comm comm = MPI_Comm_f2c(fcomm);
  if (!pname.empty()) {
    profiler = Profiler::single(pname);
  }
  instances.emplace(Instance{molpro::linalg::itsolv::create_LinearEigensystem<Rvector, Qvector, Pvector>(algorithm, options),
                             profiler, nQ, comm});
  auto& instance = instances.top();
  set_distribution(instance, range_begin, range_end);
  instance.solver->set_n_roots(nroot);
  instance.solver->set_verbosity(verbosity);
  instance.verbosity = verbosity;
  instance.has_eigenvalues = true;
  LinearEigensystemDavidson<Rvector, Qvector, Pvector>* solver =
      dynamic_cast<LinearEigensystemDavidson<Rvector, Qvector, Pvector>*>(instance.solver.get());
  if (solver) {
    solver->set_hermiticity(hermitian);
    solver->set_convergence_threshold(thresh);
    solver->set_convergence_threshold_value(thresh_value);
    //    solver_cast->propose_rspace_norm_thresh = 1.0e-14;
    //    solver_cast->set_max_size_qspace(10);
    //    solver_cast->set_reset_D(50);
    solver->logger().set_verbosity(
        verbosity > 3 ? molpro::linalg::itsolv::log::Verbosity::Trace
                      : (verbosity > 2 ? molpro::linalg::itsolv::log::Verbosity::Debug : molpro::linalg::itsolv::log::Verbosity::Info));
    solver->logger().set_min_severity(
        verbosity > 1 ? molpro::linalg::itsolv::log::Severity::Warning : molpro::linalg::itsolv::log::Severity::Error);
    solver->logger().enable_data_dumps(verbosity > 0);
  }
}

extern "C" void IterativeSolverLinearEquationsInitialize(size_t n, size_t nroot, size_t* range_begin, size_t* range_end,
                                                         const double* rhs, double aughes, double thresh,
                                                         double thresh_value, int hermitian, int verbosity,
                                                         const char* fname, int64_t fcomm, const char* algorithm,
                                                         const char* options) {
  std::shared_ptr<Profiler> profiler = nullptr;
  std::string pname(fname);
  MPI_Comm comm = MPI_Comm_f2c(fcomm);
  if (!pname.empty()) {
    profiler = Profiler::single(pname);
  }
  instances.emplace(
      Instance{molpro::linalg::itsolv::create_LinearEquations<Rvector, Qvector, Pvector>(algorithm, options), profiler,
               n, comm});
  auto& instance = instances.top();
  set_distribution(instance, range_begin, range_end);
  auto rr = CreateDistrArray(nroot, rhs);
  auto solver = dynamic_cast<LinearEquationsDavidson<Rvector, Qvector, Pvector>*>(instance.solver.get());
  if (!solver)
    throw std::runtime_error("IterativeSolverLinearEquationsInitialize: solver factory returned an unexpected type for algorithm \"" +
                             std::string(algorithm ? algorithm : "") + "\"");
  solver->set_augmented_hessian(aughes);
  solver->set_hermiticity(hermitian);
  solver->set_n_roots(nroot);
  solver->add_equations(rr);
  solver->set_convergence_threshold(thresh);
  solver->set_convergence_threshold_value(thresh_value);
  solver->logger().set_verbosity(
      verbosity > 3 ? molpro::linalg::itsolv::log::Verbosity::Trace
                    : (verbosity > 2 ? molpro::linalg::itsolv::log::Verbosity::Debug : molpro::linalg::itsolv::log::Verbosity::Info));
  solver->logger().set_min_severity(
      verbosity > 1 ? molpro::linalg::itsolv::log::Severity::Warning : molpro::linalg::itsolv::log::Severity::Error);
  solver->logger().enable_data_dumps(verbosity > 0);
  // instance.solver->m_verbosity = verbosity;
  instance.solver->set_verbosity(verbosity);
  instance.verbosity = verbosity;
}

extern "C" void IterativeSolverNonLinearEquationsInitialize(size_t n, size_t* range_begin, size_t* range_end,
                                                            double thresh, int verbosity, const char* fname,
                                                            int64_t fcomm, const char* algorithm, const char* options) {
  std::shared_ptr<Profiler> profiler = nullptr;
  std::string pname(fname);
  MPI_Comm comm = MPI_Comm_f2c(fcomm);
  if (!pname.empty()) {
    profiler = Profiler::single(pname);
  }
  instances.emplace(
      Instance{molpro::linalg::itsolv::create_NonLinearEquations<Rvector, Qvector, Pvector>(algorithm, options),
               profiler, n, comm});
  auto& instance = instances.top();
  set_distribution(instance, range_begin, range_end);
  instance.solver->set_convergence_threshold(thresh);
  // instance.solver->m_verbosity = verbosity;
  instance.solver->set_verbosity(verbosity);
  instance.verbosity = verbosity;
  molpro::linalg::itsolv::NonLinearEquationsDIIS<Rvector, Qvector, Pvector>* solver =
      dynamic_cast<molpro::linalg::itsolv::NonLinearEquationsDIIS<Rvector, Qvector, Pvector>*>(instance.solver.get());
  if (!solver)
    throw std::runtime_error(
        "IterativeSolverNonLinearEquationsInitialize: solver factory returned an unexpected type for algorithm \"" +
        std::string(algorithm ? algorithm : "") + "\"");
  solver->logger().set_verbosity(
     verbosity > 3 ? molpro::linalg::itsolv::log::Verbosity::Trace
                   : (verbosity > 2 ? molpro::linalg::itsolv::log::Verbosity::Debug : molpro::linalg::itsolv::log::Verbosity::Info));
  solver->logger().set_min_severity(
      verbosity > 1 ? molpro::linalg::itsolv::log::Severity::Warning : molpro::linalg::itsolv::log::Severity::Error);
  solver->logger().enable_data_dumps(verbosity > 0);
}

extern "C" void IterativeSolverOptimizeInitialize(size_t n, size_t* range_begin, size_t* range_end, double thresh,
                                                  double thresh_value, int verbosity, int minimize, const char* fname,
                                                  int64_t fcomm, const char* algorithm, const char* options) {
  std::shared_ptr<Profiler> profiler = nullptr;
  std::string pname(fname);
  MPI_Comm comm = MPI_Comm_f2c(fcomm);
  if (!pname.empty()) {
    profiler = Profiler::single(pname);
  }
  instances.emplace(Instance{molpro::linalg::itsolv::create_Optimize<Rvector, Qvector, Pvector>(algorithm, options),
                             profiler, n, comm});
  auto& instance = instances.top();
  set_distribution(instance, range_begin, range_end);
  instance.solver->set_n_roots(1);
  instance.solver->set_convergence_threshold(thresh);
  instance.solver->set_convergence_threshold_value(thresh_value);
  instance.solver->set_verbosity(verbosity);
  instance.verbosity = verbosity;
  molpro::linalg::itsolv::OptimizeBFGS<Rvector, Qvector, Pvector>* solver =
      dynamic_cast<molpro::linalg::itsolv::OptimizeBFGS<Rvector, Qvector, Pvector>*>(instance.solver.get());
  if (!solver)
    throw std::runtime_error(
        "IterativeSolverOptimizeInitialize: BFGS-specific configuration requested but the factory returned a different algorithm for \"" +
        std::string(algorithm ? algorithm : "") + "\"");
  solver->logger().set_verbosity(
     verbosity > 3 ? molpro::linalg::itsolv::log::Verbosity::Trace
                   : (verbosity > 2 ? molpro::linalg::itsolv::log::Verbosity::Debug : molpro::linalg::itsolv::log::Verbosity::Info));
  solver->logger().set_min_severity(
      verbosity > 1 ? molpro::linalg::itsolv::log::Severity::Warning : molpro::linalg::itsolv::log::Severity::Error);
  solver->logger().enable_data_dumps(verbosity > 0);

  instance.has_values = true;
  instance.maximise = minimize == 0;
}

extern "C" void IterativeSolverFinalize() { instances.pop(); }

extern "C" void IterativeSolverFinalizeAll() { instances.clear(); }

extern "C" int64_t IterativeSolverCurrentHandle() { return instances.current_handle(); }

extern "C" int IterativeSolverSelect(int64_t handle) { return instances.select(handle) ? 0 : 1; }

extern "C" void IterativeSolverFinalizeHandle(int64_t handle) { instances.erase(handle); }

extern "C" size_t IterativeSolverNRoots() {
  require_instance();
  return instances.top().solver->n_roots();
}

extern "C" size_t IterativeSolverDimension() {
  require_instance();
  return instances.top().dimension;
}

extern "C" void IterativeSolverAddEquation(double* rhs) {
  if (instances.empty())
    throw std::runtime_error("IterativeSolver not initialised properly");
  auto& instance = instances.top();
  if (instance.prof != nullptr)
    instance.prof->start("AddEquation");
  auto ccc = CreateDistrArray(1, rhs);
  auto* solver =
      dynamic_cast<molpro::linalg::itsolv::LinearEquationsDavidson<Rvector, Qvector, Pvector>*>(instance.solver.get());
  if (!solver)
    throw std::runtime_error("IterativeSolverAddEquation: current solver is not LinearEquationsDavidson");
  solver->add_equations(ccc[0]);
  if (instance.prof != nullptr) {
    instance.prof->stop();
  }
}

extern "C" size_t IterativeSolverAddValue(double value, double* parameters, double* action, int sync) {
  if (instances.empty())
    throw std::runtime_error("IterativeSolver not initialised properly");
  auto& instance = instances.top();
  if (instance.prof != nullptr)
    instance.prof->start("AddValue");
  auto ccc = CreateDistrArray(1, parameters);
  auto ggg = CreateDistrArray(1, action);
  if (instance.prof != nullptr)
    instance.prof->start("AddValue:Call");
  auto* solver = dynamic_cast<molpro::linalg::itsolv::Optimize<Rvector, Qvector, Pvector>*>(instance.solver.get());
  if (!solver)
    throw std::runtime_error("IterativeSolverAddValue: current solver is not an Optimize solver");
  if (instance.maximise) {
    value = -value;
    ggg[0].scal(-1);
  }
  size_t working_set_size = solver->add_vector(ccc[0], ggg[0], value) > 0 ? 1 : 0;
  // the caller preconditions the returned residual with its own (negative) Hessian diagonal, which gives the step
  // that minimising the negative objective needs: (-g)/(-d) = g/d
  if (instance.maximise)
    ggg[0].scal(-1);
  if (instance.prof != nullptr) {
    instance.prof->stop();
    instance.prof->start("AddValue:Sync");
  }
  if (sync) {
    DistrArraySynchronize(1, ccc, parameters);
    DistrArraySynchronize(1, ggg, action);
  }
  if (instance.prof != nullptr) {
    instance.prof->stop();
    instance.prof->stop();
  }
  return working_set_size;
}

extern "C" size_t IterativeSolverAddVector(size_t buffer_size, double* parameters, double* action, int sync) {
  if (instances.empty())
    throw std::runtime_error("IterativeSolver not initialised properly");
  auto& instance = instances.top();
  if (instance.prof != nullptr)
    instance.prof->start("AddVector");
  auto cc = CreateDistrArray(buffer_size, parameters);
  auto gg = CreateDistrArray(buffer_size, action);
  if (instance.prof != nullptr)
    instance.prof->start("AddVector:Call");
  auto working_set_size = instance.solver->add_vector(cc, gg);
  if (instance.prof != nullptr) {
    instance.prof->stop();
    instance.prof->start("AddVector:Sync");
  }
  if (sync) {
    DistrArraySynchronize(instance.solver->working_set().size(), cc, parameters);
    DistrArraySynchronize(instance.solver->working_set().size(), gg, action);
  }
  if (instance.prof != nullptr) {
    instance.prof->stop();
    instance.prof->stop();
  }
  return working_set_size;
}

extern "C" void IterativeSolverSolution(int nroot, int* roots, double* parameters, double* action, int sync) {
  if (instances.empty())
    throw std::runtime_error("IterativeSolver not initialised properly");
  auto& instance = instances.top();
  if (instance.prof != nullptr)
    instance.prof->start("Solution");
  auto cc = CreateDistrArray(nroot, parameters);
  auto gg = CreateDistrArray(nroot, action);
  std::vector<int> croots;
  for (int i = 0; i < nroot; i++) {
    croots.push_back(*(roots + i));
  }
  if (instance.prof != nullptr)
    instance.prof->start("Solution:Call");
  instance.solver->solution(croots, cc, gg);
  if (instance.maximise)
    for (auto& g : gg)
      g.scal(-1);
  if (instance.prof != nullptr) {
    instance.prof->stop();
    instance.prof->start("Solution:Sync");
  }
  if (sync) {
    DistrArraySynchronize(nroot, cc, parameters);
    DistrArraySynchronize(nroot, gg, action);
  }
  if (instance.prof != nullptr) {
    instance.prof->stop();
    instance.prof->stop();
  }
}

extern "C" size_t IterativeSolverEndIteration(size_t buffer_size, double* solution, double* residual, int sync) {
  if (instances.empty())
    throw std::runtime_error("IterativeSolver not initialised properly");
  auto& instance = instances.top();
  if (instance.prof != nullptr)
    instance.prof->start("EndIter");
  auto cc = CreateDistrArray(buffer_size, solution);
  auto gg = CreateDistrArray(buffer_size, residual);
  if (instance.prof != nullptr)
    instance.prof->start("EndIter:Call");
  auto result = instance.solver->end_iteration(cc, gg);
  if (instance.prof != nullptr) {
    instance.prof->stop();
    instance.prof->start("EndIter:Sync");
  }
  if (sync) {
    DistrArraySynchronize(instance.solver->working_set().size(), cc, solution);
    DistrArraySynchronize(instance.solver->working_set().size(), gg, residual);
  }
  if (instance.prof != nullptr) {
    instance.prof->stop();
    instance.prof->stop();
  }
  return result;
}

extern "C" int IterativeSolverEndIterationNeeded(){
  if (instances.empty())
    throw std::runtime_error("IterativeSolver not initialised properly");
  auto& instance = instances.top();
  return instance.solver->end_iteration_needed() ? 1:0;
}

extern "C" size_t IterativeSolverAddP(size_t buffer_size, size_t nP, const size_t* offsets, const size_t* indices,
                                      const double* coefficients, const double* pp, double* parameters, double* action,
                                      int sync, apply_on_p_t func) {
  if (instances.empty())
    throw std::runtime_error("IterativeSolver not initialised properly");
  auto& instance = instances.top();
  instance.apply_on_p_fort = func;
  if (instance.prof != nullptr)
    instance.prof->start("AddP");
  auto cc = CreateDistrArray(buffer_size, parameters);
  auto gg = CreateDistrArray(buffer_size, action);
  std::vector<Pvector> Pvectors;
  Pvectors.reserve(nP);
  for (size_t p = 0; p < nP; p++) {
    Pvector ppp;
    for (size_t k = offsets[p]; k < offsets[p + 1]; k++)
      ppp.insert(std::pair<size_t, Rvector::value_type>(indices[k], coefficients[k]));
    Pvectors.emplace_back(ppp);
  }
  using vectorP = std::vector<double>;
  using molpro::linalg::itsolv::CVecRef;
  using molpro::linalg::itsolv::VecRef;
  std::function<void(const std::vector<vectorP>&, const CVecRef<Pvector>&, const VecRef<Rvector>&)> apply_on_p =
      apply_on_p_c;
  if (instance.prof != nullptr)
    instance.prof->start("AddP:Call");
  size_t working_set_size = instance.solver->add_p(
      cwrap(Pvectors),
      Span<Rvector::value_type>(&const_cast<double*>(pp)[0], (instance.solver->dimensions().oP() + nP) * nP), wrap(cc),
      wrap(gg), apply_on_p);
  if (instance.prof != nullptr) {
    instance.prof->stop();
    instance.prof->start("AddP:Sync");
  }
  if (sync) {
    DistrArraySynchronize(working_set_size, cc, parameters);
    DistrArraySynchronize(working_set_size, gg, action);
  }
  if (instance.prof != nullptr) {
    instance.prof->stop();
    instance.prof->stop();
  }
  return working_set_size;
}

namespace {
void require_instance() {
  if (instances.empty())
    throw std::runtime_error("IterativeSolver not initialised properly");
}
} // namespace

extern "C" void IterativeSolverErrors(double* errors) {
  require_instance();
  auto& instance = instances.top();
  size_t k = 0;
  for (const auto& e : instance.solver.get()->errors())
    errors[k++] = e;
  return;
}

extern "C" void IterativeSolverEigenvalues(double* eigenvalues) {
  require_instance();
  auto& instance = instances.top();
  size_t k = 0;
  LinearEigensystem<Rvector, Qvector, Pvector>* solver_cast =
      dynamic_cast<LinearEigensystem<Rvector, Qvector, Pvector>*>(instance.solver.get());
  if (solver_cast) {
    for (const auto& e : solver_cast->eigenvalues())
      eigenvalues[k++] = e;
  }
}

extern "C" void IterativeSolverWorkingSetEigenvalues(double* eigenvalues) {
  require_instance();
  auto& instance = instances.top();
  size_t k = 0;
  LinearEigensystemDavidson<Rvector, Qvector, Pvector>* solver_cast =
      dynamic_cast<LinearEigensystemDavidson<Rvector, Qvector, Pvector>*>(instance.solver.get());
  if (solver_cast) {
    for (const auto& e : solver_cast->working_set_eigenvalues())
      eigenvalues[k++] = e;
  }
}

extern "C" size_t IterativeSolverSuggestP(const double* solution, const double* residual, size_t maximumNumber,
                                          double threshold, size_t* indices) {
  require_instance();
  auto& instance = instances.top();
  if (instance.prof != nullptr)
    instance.prof->start("SuggestP");
  auto cc = CreateDistrArray(instance.solver->n_roots(), solution);
  auto gg = CreateDistrArray(instance.solver->n_roots(), residual);
  auto result = instance.solver->suggest_p(cwrap(cc), molpro::linalg::itsolv::cwrap(gg), maximumNumber, threshold);
  for (size_t i = 0; i < result.size(); i++) {
    indices[i] = result[i];
  }
  if (instance.prof != nullptr)
    instance.prof->stop();
  return result.size();
}

extern "C" void IterativeSolverPrintStatistics() {
  require_instance();
  molpro::cout << instances.top().solver->statistics() << std::endl;
}

extern "C" int IterativeSolverConverged() {
  require_instance();
  const auto& solver = instances.top().solver;
  // as in solve(): an empty working set is not enough, the errors must also be within the threshold
  const auto& errors = solver->errors();
  return solver->working_set().empty() and not errors.empty() and
                 *std::max_element(errors.begin(), errors.end()) <= solver->convergence_threshold()
             ? 1
             : 0;
}

extern "C" int IterativeSolverNonLinear() {
  require_instance();
  return instances.top().solver->nonlinear() ? 1 : 0;
}
extern "C" int IterativeSolverHasValues() {
  require_instance();
  return instances.top().has_values ? 1 : 0;
}
extern "C" int IterativeSolverHasEigenvalues() {
  require_instance();
  return instances.top().has_eigenvalues ? 1 : 0;
}

extern "C" void IterativeSolverSetDiagonals(const double* diagonals) {
  require_instance();
  instances.top().diagonals.reset(new Qvector(CreateDistrArray(1, diagonals).front()));
}
extern "C" void IterativeSolverDiagonals(double* diagonals) {
  require_instance();
  CreateDistrArray(1, diagonals).front().copy(*instances.top().diagonals);
}
extern "C" double IterativeSolverValue() {
  require_instance();
  return instances.top().maximise ? -instances.top().solver->value() : instances.top().solver->value();
}
// The initialisers set the logger with their own mapping from the requested print level, which cannot be inverted
// (levels 0 to 2 all give log::Verbosity::Info), so return the level that was requested.
extern "C" int IterativeSolverVerbosity() {
  require_instance();
  return instances.top().verbosity;
}
extern "C" int IterativeSolverMaxIter() {
  require_instance();
  return instances.top().solver->get_max_iter();
}
extern "C" void IterativeSolverSetMaxIter(int max_iter) {
  require_instance();
  instances.top().solver->set_max_iter(max_iter);
}
/*!
 * @brief C binding of mpi::comm_global(), suitable for calling from Fortran
 */
extern "C" int64_t IterativeSolver_mpicomm_global() { return (int64_t)MPI_Comm_c2f(molpro::mpi::comm_global()); }

/*!
 * @brief C binding of mpi::comm_self(), suitable for calling from Fortran
 */
extern "C" int64_t IterativeSolver_mpicomm_self() { return (int64_t)MPI_Comm_c2f(molpro::mpi::comm_self()); }

/*!
 * @brief C binding of mpi::size_global(), suitable for calling from Fortran
 */
extern "C" int64_t IterativeSolver_mpi_size_global() { return molpro::mpi::size_global(); }

/*!
 * @brief C binding of mpi::rank_global(), suitable for calling from Fortran
 */
extern "C" int64_t IterativeSolver_mpi_rank_global() { return molpro::mpi::rank_global(); }

/*!
 * @brief C binding of mpi::init(), suitable for calling from Fortran
 */
extern "C" int IterativeSolver_mpi_init() { return molpro::mpi::init(); }

/*!
 * @brief C binding of mpi::finalize(), suitable for calling from Fortran
 */
extern "C" int IterativeSolver_mpi_finalize() { return molpro::mpi::finalize(); }
