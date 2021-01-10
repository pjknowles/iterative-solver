#include "IterativeSolverC.h"
#include "molpro/ProfilerSingle.h"
#include <memory>
#include <mpi.h>
#include <stack>
#include <string>
#include <tuple>
#ifdef HAVE_PPIDD_H
#include <ppidd.h>
#endif
#ifdef LINEARALGEBRA_ARRAY_GA
#include "ga-mpi.h"
#include "ga.h"
#endif

#ifdef LINEARALGEBRA_ARRAY_HDF5
#include <molpro/linalg/array/DistrArrayHDF5.h>
#else
#include <molpro/linalg/array/DistrArrayFile.h>
#endif
#include <molpro/linalg/array/DistrArrayMPI3.h>
#include <molpro/linalg/array/Span.h>
#include <molpro/linalg/array/util/Distribution.h>
#include <molpro/linalg/array/util/gather_all.h>
#include <molpro/linalg/itsolv/ArrayHandlers.h>
#include <molpro/linalg/itsolv/LinearEigensystemDavidson.h>
#include <molpro/linalg/itsolv/LinearEquationsDavidson.h>
#include <molpro/linalg/itsolv/SolverFactory.h>

using molpro::Profiler;
using molpro::linalg::array::Span;
using molpro::linalg::array::util::gather_all;
using molpro::linalg::itsolv::ArrayHandlers;
using molpro::linalg::itsolv::LinearEigensystem;
using molpro::linalg::itsolv::LinearEquations;
using molpro::linalg::itsolv::NonLinearEquations;
using molpro::linalg::itsolv::Optimize;
using molpro::linalg::itsolv::IterativeSolver;
using molpro::linalg::itsolv::LinearEigensystemDavidson;
using molpro::linalg::itsolv::LinearEquationsDavidson;

using Rvector = molpro::linalg::array::DistrArrayMPI3;
#ifdef LINEARALGEBRA_ARRAY_HDF5
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

// FIXME Only top solver is active at any one time. This should be documented somewhere.

namespace {
typedef void (*Apply_on_p_fort)(const double*, double*, const size_t, const size_t*);
struct Instance {
  Instance(std::unique_ptr<IterativeSolver<Rvector, Qvector, Pvector>> solver, std::shared_ptr<Profiler> prof,
           size_t dimension, MPI_Comm comm)
      : solver(std::move(solver)), prof(std::move(prof)), dimension(dimension), comm(comm){};
  std::unique_ptr<IterativeSolver<Rvector, Qvector, Pvector>> solver;
  std::shared_ptr<Profiler> prof;
  Apply_on_p_fort apply_on_p_fort;
  size_t dimension;
  MPI_Comm comm;
};
std::stack<Instance> instances;
} // namespace

std::vector<Rvector> CreateDistrArray(size_t nvec, double* data) {
  auto& instance = instances.top();
  MPI_Comm ccomm = instance.comm;
  int mpi_rank;
  MPI_Comm_rank(ccomm, &mpi_rank);
  std::vector<Rvector> c;
  c.reserve(nvec);
  for (size_t ivec = 0; ivec < nvec; ivec++) {
    c.emplace_back(instance.dimension, ccomm);
    auto crange = c.back().distribution().range(mpi_rank);
    auto clength = crange.second - crange.first;
    c.back().allocate_buffer(
        Span<typename Rvector::value_type>(&data[ivec * instance.dimension + crange.first], clength));
  }
  return c;
}

std::vector<Rvector> CreateDistrArrayConst(size_t nvec, const double* data) {
  auto& instance = instances.top();
  MPI_Comm ccomm = instance.comm;
  int mpi_rank;
  MPI_Comm_rank(ccomm, &mpi_rank);
  std::vector<Rvector> c;
  c.reserve(nvec);
  for (size_t ivec = 0; ivec < nvec; ivec++) {
    c.emplace_back(instance.dimension, ccomm);
    auto crange = c.back().distribution().range(mpi_rank);
    auto clength = crange.second - crange.first;
    c.back().allocate_buffer(
        Span<typename Rvector::value_type>(&const_cast<double*>(data)[ivec * instance.dimension + crange.first], clength));
  }
  return c;
}

void synchronize(size_t nvec, std::vector<Rvector>& c, double* data) {
  auto& instance = instances.top();
  MPI_Comm ccomm = instance.comm;
  for (size_t ivec = 0; ivec < nvec; ivec++) {
    gather_all(c[ivec].distribution(), ccomm, &data[ivec * instance.dimension]);
  }
}

void apply_on_p_c(const std::vector<vectorP>& pvectors, const CVecRef<Pvector>& pspace, const VecRef<Rvector>& action) {
  auto& instance = instances.top();
  MPI_Comm ccomm = instance.comm;
  int mpi_rank, mpi_size;
  MPI_Comm_rank(ccomm, &mpi_rank);
  MPI_Comm_size(ccomm, &mpi_size);
  std::vector<size_t> ranges;
  size_t update_size = pvectors.size();
  ranges.reserve(update_size * 2);
  for (size_t k = 0; k < update_size; ++k) {
    auto range = action[k].get().distribution().range(mpi_rank);
    ranges.push_back(range.first);
    ranges.push_back(range.second);
  }
  std::vector<double> pvecs_to_send;
  for (size_t i = 0; i < update_size; i++) {
    for (auto j : pvectors[i]) {
      pvecs_to_send.push_back(j);
    }
  }
  instance.apply_on_p_fort(pvecs_to_send.data(), &(*action.front().get().local_buffer())[0], update_size,
                           ranges.data());
}

extern "C" void IterativeSolverLinearEigensystemInitialize(size_t nQ, size_t nroot, size_t* range_begin,
                                                           size_t* range_end, double thresh, double thresh_value,
                                                           int hermitian, int verbosity, const char* fname,
                                                           int64_t fcomm, const char* algorithm) {
  std::shared_ptr<Profiler> profiler = nullptr;
  std::string pname(fname);
  MPI_Comm comm = MPI_Comm_f2c(fcomm);
  if (!pname.empty()) {
    profiler = molpro::ProfilerSingle::instance(pname, comm);
  }
  auto handlers = std::make_shared<ArrayHandlers<Rvector, Qvector, Pvector>>();
  instances.emplace(Instance{molpro::linalg::itsolv::create_LinearEigensystem<Rvector, Qvector, Pvector>(algorithm, "")
                             //              std::make_unique<LinearEigensystem<Rvector, Qvector, Pvector>>(handlers)
                             ,
                             profiler, nQ, comm});
  auto& instance = instances.top();
  instance.solver->set_n_roots(nroot);
  LinearEigensystemDavidson<Rvector, Qvector, Pvector>* solver_cast =
      dynamic_cast<LinearEigensystemDavidson<Rvector, Qvector, Pvector>*>(instance.solver.get());
  if (solver_cast) {
    solver_cast->set_hermiticity(hermitian);
    solver_cast->set_convergence_threshold(thresh);
    solver_cast->set_convergence_threshold_value(thresh_value);
    //    solver_cast->propose_rspace_norm_thresh = 1.0e-14;
    //    solver_cast->set_max_size_qspace(10);
    //    solver_cast->set_reset_D(50);
    solver_cast->logger->max_trace_level = molpro::linalg::itsolv::Logger::None;
    solver_cast->logger->max_warn_level = molpro::linalg::itsolv::Logger::Error;
    solver_cast->logger->data_dump = false;
  }
  //TODO: we just need range here, no reason to allocate vectors
  std::vector<Rvector> x;
  std::vector<Rvector> g;
  for (size_t root = 0; root < instance.solver->n_roots(); root++) {
    x.emplace_back(nQ, comm);
    g.emplace_back(nQ, comm);
  }
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);
  std::tie(*range_begin, *range_end) = x[0].distribution().range(mpi_rank);
}

extern "C" void IterativeSolverLinearEquationsInitialize(size_t n, size_t nroot, size_t* range_begin, size_t* range_end,
                                                         const double* rhs, double aughes, double thresh,
                                                         double thresh_value, int hermitian, int verbosity,
                                                         const char* fname, int64_t fcomm, const char* algorithm) {
  std::shared_ptr<Profiler> profiler = nullptr;
  std::string pname(fname);
  MPI_Comm comm = MPI_Comm_f2c(fcomm);
  if (!pname.empty()) {
    profiler = molpro::ProfilerSingle::instance(pname, comm);
  }
  auto handlers = std::make_shared<ArrayHandlers<Rvector, Qvector, Pvector>>();
  instances.emplace(
      Instance{molpro::linalg::itsolv::create_LinearEquations<Rvector, Qvector, Pvector>(algorithm, "")
               //        std::make_unique<LinearEquations<Rvector, Qvector, Pvector>>(rr, handlers, aughes)
               ,
               profiler, n, comm});
  auto& instance = instances.top();
  auto rr = CreateDistrArrayConst(nroot, rhs);
  auto solver = dynamic_cast<LinearEquations<Rvector, Qvector, Pvector>*>(instance.solver.get());
  auto solverLinearEquations = dynamic_cast<LinearEquationsDavidson<Rvector, Qvector, Pvector>*>(instance.solver.get());
  solver->set_n_roots(nroot);
  solverLinearEquations->add_equations(rr);
  solver->set_convergence_threshold(thresh);
  solver->set_convergence_threshold_value(thresh_value);
  // instance.solver->m_verbosity = verbosity;
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);
  std::tie(*range_begin, *range_end) = rr[0].distribution().range(mpi_rank);
}

extern "C" void IterativeSolverNonLinearEquationsInitialize(size_t n, size_t* range_begin, size_t* range_end,
                                                            double thresh, int verbosity, const char* fname,
                                                            int64_t fcomm, const char* algorithm) {
  std::shared_ptr<Profiler> profiler = nullptr;
  std::string pname(fname);
  MPI_Comm comm = MPI_Comm_f2c(fcomm);
  if (!pname.empty()) {
    profiler = molpro::ProfilerSingle::instance(pname, comm);
  }
  auto handlers = std::make_shared<ArrayHandlers<Rvector, Qvector, Pvector>>();
  instances.emplace(Instance{
      molpro::linalg::itsolv::create_NonLinearEquations<Rvector, Qvector, Pvector>(algorithm, ""), profiler, n, comm});
  auto& instance = instances.top();
  instance.solver->set_convergence_threshold(thresh);
  // instance.solver->m_verbosity = verbosity;
  Rvector x(n, comm);
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);
  std::tie(*range_begin, *range_end) = x.distribution().range(mpi_rank);
}

extern "C" void IterativeSolverOptimizeInitialize(size_t n, size_t* range_begin, size_t* range_end, double thresh,
                                                  double thresh_value, int verbosity, int minimize, const char* fname,
                                                  int64_t fcomm, const char* algorithm) {
  std::shared_ptr<Profiler> profiler = nullptr;
  std::string pname(fname);
  MPI_Comm comm = MPI_Comm_f2c(fcomm);
  if (!pname.empty()) {
    profiler = molpro::ProfilerSingle::instance(pname, comm);
  }
  auto handlers = std::make_shared<ArrayHandlers<Rvector, Qvector, Pvector>>();
  instances.emplace(
      Instance{molpro::linalg::itsolv::create_Optimize<Rvector, Qvector, Pvector>(algorithm, ""), profiler, n, comm});
  auto& instance = instances.top();
  instance.solver->set_n_roots(1);
  instance.solver->set_convergence_threshold(thresh);
  instance.solver->set_convergence_threshold_value(thresh_value);
  // instance.solver->m_verbosity = verbosity;
  Rvector x(n, comm);
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);
  std::tie(*range_begin, *range_end) = x.distribution().range(mpi_rank);
}

extern "C" void IterativeSolverFinalize() { instances.pop(); }

extern "C" size_t IterativeSolverAddValue(double value, double* parameters, double* action, int sync) {
  auto& instance = instances.top();
  auto ccc = CreateDistrArray(1, parameters);
  auto ggg = CreateDistrArray(1, action);
  size_t working_set_size =
      dynamic_cast<molpro::linalg::itsolv::Optimize<Rvector, Qvector, Pvector>*>(instance.solver.get())
              ->add_value(ccc[0], value, ggg[0])
          ? 1
          : 0;
  if (sync) {
    synchronize(1, ccc, parameters);
    synchronize(1, ggg, action);
  }
  return working_set_size;
}

extern "C" size_t IterativeSolverAddVector(size_t buffer_size, double* parameters, double* action, int sync) {
  if (instances.empty())
    throw std::runtime_error("IterativeSolver not initialised properly");
  auto& instance = instances.top();
  if (instance.prof != nullptr)
    instance.prof->start("AddVector");
  size_t working_set_size = instance.solver->working_set().size();
  auto cc = CreateDistrArray(working_set_size, parameters);
  auto gg = CreateDistrArray(working_set_size, action);
  if (instance.prof != nullptr)
    instance.prof->start("AddVector:Update");
  working_set_size = instance.solver->add_vector(cc, gg);
  if (instance.prof != nullptr)
    instance.prof->stop("AddVector:Update");

  if (instance.prof != nullptr)
    instance.prof->start("AddVector:Sync");
  if (sync) {
    synchronize(instance.solver->working_set().size(), cc, parameters);
    synchronize(instance.solver->working_set().size(), gg, action);
  }
  if (instance.prof != nullptr)
    instance.prof->stop("AddVector:Sync");
  if (instance.prof != nullptr)
    instance.prof->stop("AddVector");
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
  if (instance.prof != nullptr)
    instance.prof->stop("Solution:Call");

  if (instance.prof != nullptr)
    instance.prof->start("Solution:Sync");
  if (sync) {
    synchronize(nroot, cc, parameters);
    synchronize(nroot, gg, action);
  }
  if (instance.prof != nullptr)
    instance.prof->stop("Solution:Sync");
  if (instance.prof != nullptr)
    instance.prof->stop("Solution");
}

extern "C" int IterativeSolverEndIteration(size_t buffer_size, double* solution, double* residual, int sync) {
  if (instances.empty())
    throw std::runtime_error("IterativeSolver not initialised properly");
  auto& instance = instances.top();
  if (instance.prof != nullptr)
    instance.prof->start("EndIter");
  auto cc = CreateDistrArray(buffer_size, solution);
  auto gg = CreateDistrArray(buffer_size, residual);
  if (instance.prof != nullptr)
    instance.prof->start("EndIter:Call");
  int result = instance.solver->end_iteration(cc, gg);
  if (instance.prof != nullptr)
    instance.prof->stop("EndIter:Call");
  if (instance.prof != nullptr)
    instance.prof->start("AddVector:Sync");
  if (sync) { // should be over instance.solver->working_set().size()
    synchronize(instance.solver->working_set().size(), cc, solution);
    synchronize(instance.solver->working_set().size(), gg, residual);
  }
  if (instance.prof != nullptr)
    instance.prof->stop("AddVector:Sync");
  if (instance.prof != nullptr)
    instance.prof->stop("EndIter");
  return result;
}

extern "C" size_t IterativeSolverAddP(size_t buffer_size, size_t nP, const size_t* offsets, const size_t* indices,
                                      const double* coefficients, const double* pp, double* parameters, double* action,
                                      int sync, void (*func)(const double*, double*, const size_t, const size_t*)) {
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
      molpro::linalg::itsolv::cwrap(Pvectors),
      Span<Rvector::value_type>(&const_cast<double*>(pp)[0], (instance.solver->dimensions().oP + nP) * nP),
      molpro::linalg::itsolv::wrap(cc), molpro::linalg::itsolv::wrap(gg), apply_on_p);
  if (instance.prof != nullptr)
    instance.prof->stop("AddP:Call");
  if (instance.prof != nullptr)
    instance.prof->start("AddP:Sync");
  if (sync) {
    synchronize(working_set_size, cc, parameters);
    synchronize(working_set_size, gg, action);
  }
  if (instance.prof != nullptr)
    instance.prof->stop("AddP:Sync");
  if (instance.prof != nullptr)
    instance.prof->stop("AddP");
  return working_set_size;
}

extern "C" void IterativeSolverErrors(double* errors) {
  auto& instance = instances.top();
  size_t k = 0;
  for (const auto& e : instance.solver.get()->errors())
    errors[k++] = e;
  return;
}

extern "C" void IterativeSolverEigenvalues(double* eigenvalues) {
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
  auto& instance = instances.top();
  if (instance.prof != nullptr)
    instance.prof->start("SuggestP");
  auto cc = CreateDistrArrayConst(instance.solver->n_roots(), solution);
  auto gg = CreateDistrArrayConst(instance.solver->n_roots(), residual);
  auto result = instance.solver->suggest_p(molpro::linalg::itsolv::cwrap(cc), molpro::linalg::itsolv::cwrap(gg),
                                           maximumNumber, threshold);
  for (size_t i = 0; i < result.size(); i++) {
    indices[i] = result[i];
  }
  if (instance.prof != nullptr)
    instance.prof->stop("SuggestP");
  return result.size();
}

extern "C" void IterativeSolverPrintStatistics() { molpro::cout << instances.top().solver->statistics() << std::endl; }

extern "C" int64_t mpicomm_self() {
  int flag;
  MPI_Initialized(&flag);
  if (!flag)
    return 0;
  return MPI_Comm_c2f(MPI_COMM_SELF);
}

extern "C" int64_t mpicomm_global() {
  int flag;
  MPI_Initialized(&flag);
  if (!flag) {
    MPI_Init(0, nullptr);
    return MPI_Comm_c2f(MPI_COMM_WORLD);
  }
#ifdef HAVE_PPIDD_H
  {
    int64_t size;
    PPIDD_Size(&size);
    if (size > 0)
      return PPIDD_Worker_comm();
  }
#else
#ifdef LINEARALGEBRA_ARRAY_GA
  if (GA_MPI_Comm() != NULL && GA_MPI_Comm() != MPI_COMM_NULL) {
    return MPI_Comm_c2f(GA_MPI_Comm());
  }
#endif
#endif
  return MPI_Comm_c2f(MPI_COMM_WORLD);
}
