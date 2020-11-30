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

#include <molpro/linalg/array/DistrArrayHDF5.h>
#include <molpro/linalg/array/DistrArrayMPI3.h>
#include <molpro/linalg/array/Span.h>
#include <molpro/linalg/array/util/Distribution.h>
#include <molpro/linalg/array/util/gather_all.h>
#include <molpro/linalg/itsolv/ArrayHandlers.h>
#include <molpro/linalg/itsolv/LinearEigensystem.h>

using molpro::Profiler;
using molpro::linalg::array::Span;
using molpro::linalg::array::util::gather_all;
using molpro::linalg::itsolv::ArrayHandlers;
using molpro::linalg::itsolv::ILinearEigensystem;
using molpro::linalg::itsolv::ILinearEquations;
using molpro::linalg::itsolv::INonLinearEquations;
using molpro::linalg::itsolv::IOptimize;
using molpro::linalg::itsolv::IterativeSolver;
using molpro::linalg::itsolv::LinearEigensystem;

using Rvector = molpro::linalg::array::DistrArrayMPI3;
using Qvector = molpro::linalg::array::DistrArrayMPI3;
using Pvector = std::map<size_t, double>;

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

extern "C" void IterativeSolverLinearEigensystemInitialize(size_t n, size_t nroot, size_t range_begin, size_t range_end,
                                                           double thresh, int hermitian, int verbosity,
                                                           const char* fname, int64_t fcomm, int lmppx) {
  std::shared_ptr<Profiler> profiler = nullptr;
  std::string pname(fname);
  int flag;
  MPI_Initialized(&flag);
  MPI_Comm comm;
  if (!flag) {
#ifdef HAVE_PPIDD_H
    PPIDD_Initialize(0, nullptr, PPIDD_IMPL_DEFAULT);
    comm = MPI_Comm_f2c(PPIDD_Worker_comm());
#else
    MPI_Init(0, nullptr);
    comm = MPI_COMM_WORLD;
#endif
  } else if (lmppx != 0) {
    comm = MPI_COMM_SELF;
  } else {
    // TODO: Check this is safe. Will crash if handle is invalid.
    comm = MPI_Comm_f2c(fcomm);
  }
  // TODO: what if lmppx != 0 ?
  if (!pname.empty()) {
    profiler = molpro::ProfilerSingle::instance(pname, comm);
  }
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);
  auto handlers = std::make_shared<ArrayHandlers<Rvector, Qvector, Pvector>>();
  instances.emplace(
      Instance{std::make_unique<LinearEigensystem<Rvector, Qvector, Pvector>>(handlers), profiler, n, comm});
  auto& instance = instances.top();
  instance.solver->set_n_roots(nroot);
  instance.solver->set_hermitian(hermitian);
  LinearEigensystem<Rvector, Qvector, Pvector>* solver_cast =
      dynamic_cast<LinearEigensystem<Rvector, Qvector, Pvector>*>(instance.solver.get());
  if (solver_cast) {
    solver_cast->set_convergence_threshold(1.0e-12);
    solver_cast->propose_rspace_norm_thresh = 1.0e-14;
    solver_cast->set_max_size_qspace(10);
    solver_cast->set_reset_D(50);
    solver_cast->logger->max_trace_level = molpro::linalg::itsolv::Logger::None;
    solver_cast->logger->max_warn_level = molpro::linalg::itsolv::Logger::Error;
    solver_cast->logger->data_dump = false;
  }
  std::vector<Rvector> x;
  std::vector<Rvector> g;
  for (size_t root = 0; root < instance.solver->n_roots(); root++) {
    x.emplace_back(n, comm);
    g.emplace_back(n, comm);
  }
  std::tie(range_begin, range_end) = x[0].distribution().range(mpi_rank);
}

extern "C" void IterativeSolverLinearEquationsInitialize(size_t n, size_t nroot, size_t range_begin, size_t range_end,
                                                         const double* rhs, double aughes, double thresh, int hermitian,
                                                         int verbosity, const char* fname, int64_t fcomm, int lmppx) {
  /*
    std::shared_ptr<Profiler> profiler = nullptr;
    std::string pname(fname);
    int flag;
    MPI_Initialized(&flag);
    MPI_Comm comm;
    if (!flag) {
  #ifdef HAVE_PPIDD_H
      PPIDD_Initialize(0, nullptr, PPIDD_IMPL_DEFAULT);
      comm = MPI_Comm_f2c(PPIDD_Worker_comm());
  #else
      MPI_Init(0, nullptr);
      comm = MPI_COMM_WORLD;
  #endif
    } else if (lmppx != 0) {
      comm = MPI_COMM_SELF;
    } else {
      comm = MPI_Comm_f2c(fcomm);
    }
    if (!pname.empty()) {
      profiler = molpro::ProfilerSingle::instance(pname, comm);
    }
    int mpi_rank;
    MPI_Comm_rank(comm, &mpi_rank);
    auto handlers = std::make_shared<ArrayHandlers<Rvector, Qvector, Pvector>>();
    std::vector<Rvector> rr;
    rr.reserve(nroot);
    for (size_t root = 0; root < nroot; root++) {
      rr.emplace_back(n, comm);
      auto rrrange = rr.back().distribution().range(mpi_rank);
      auto rrn = rrrange.second - rrrange.first;
      rr.back().allocate_buffer(Span<Rvector::value_type>(&const_cast<double*>(rhs)[root * n + rrrange.first], rrn));
    }
    instances.emplace(
        Instance{std::make_unique<LinearEquations<Rvector, Qvector, Pvector>>(rr, handlers, aughes), profiler, n,
  comm}); auto& instance = instances.top(); instance.solver->set_n_roots(nroot);
    //instance.solver->m_thresh = thresh;
    //instance.solver->m_verbosity = verbosity;
    std::tie(range_begin, range_end) = rr[0].distribution().range(mpi_rank);
    */
}

extern "C" void IterativeSolverDIISInitialize(size_t n, size_t range_begin, size_t range_end, double thresh,
                                              int verbosity, const char* fname, int64_t fcomm, int lmppx) {
  /*
    std::shared_ptr<Profiler> profiler = nullptr;
    std::string pname(fname);
    int flag;
    MPI_Initialized(&flag);
    MPI_Comm comm;
    if (!flag) {
  #ifdef HAVE_PPIDD_H
      PPIDD_Initialize(0, nullptr, PPIDD_IMPL_DEFAULT);
      comm = MPI_Comm_f2c(PPIDD_Worker_comm());
  #else
      MPI_Init(0, nullptr);
      comm = MPI_COMM_WORLD;
  #endif
    } else if (lmppx != 0) {
      comm = MPI_COMM_SELF;
    } else {
      comm = MPI_Comm_f2c(fcomm);
    }
    if (!pname.empty()) {
      profiler = molpro::ProfilerSingle::instance(pname, comm);
    }
    int mpi_rank;
    MPI_Comm_rank(comm, &mpi_rank);
    auto handlers = std::make_shared<ArrayHandlers<Rvector, Qvector, Pvector>>();
    instances.emplace(Instance{std::make_unique<DIIS<Rvector, Qvector, Pvector>>(handlers), profiler, n, comm});
    auto& instance = instances.top();
    //instance.solver->m_thresh = thresh;
    //instance.solver->m_verbosity = verbosity;
    Rvector x(n, comm);
    std::tie(range_begin, range_end) = x.distribution().range(mpi_rank);
    */
}

extern "C" void IterativeSolverOptimizeInitialize(size_t n, size_t range_begin, size_t range_end, double thresh,
                                                  int verbosity, char* algorithm, int minimize, const char* fname,
                                                  int64_t fcomm, int lmppx) {
  /*
    std::shared_ptr<Profiler> profiler = nullptr;
    std::string pname(fname);
    int flag;
    MPI_Initialized(&flag);
    MPI_Comm comm;
    if (!flag) {
  #ifdef HAVE_PPIDD_H
      PPIDD_Initialize(0, nullptr, PPIDD_IMPL_DEFAULT);
      comm = MPI_Comm_f2c(PPIDD_Worker_comm());
  #else
      MPI_Init(0, nullptr);
      comm = MPI_COMM_WORLD;
  #endif
    } else if (lmppx != 0) {
      comm = MPI_COMM_SELF;
    } else {
      comm = MPI_Comm_f2c(fcomm);
    }
    if (!pname.empty()) {
      profiler = molpro::ProfilerSingle::instance(pname, comm);
    }
    int mpi_rank;
    MPI_Comm_rank(comm, &mpi_rank);
    auto handlers = std::make_shared<ArrayHandlers<Rvector, Qvector, Pvector>>();
    if (*algorithm)
      instances.emplace(
          Instance{std::make_unique<Optimize<Rvector, Qvector, Pvector>>(handlers, algorithm, minimize != 0), profiler,
  n, comm}); else instances.emplace(Instance{std::make_unique<Optimize<Rvector, Qvector, Pvector>>(handlers), profiler,
  n, comm}); auto& instance = instances.top(); instance.solver->set_n_roots(1);
    //instance.solver->m_thresh = thresh;
    //instance.solver->m_verbosity = verbosity;
    Rvector x(n, comm);
    std::tie(range_begin, range_end) = x.distribution().range(mpi_rank);
    */
}

extern "C" void IterativeSolverFinalize() { instances.pop(); }

extern "C" size_t IterativeSolverAddValue(double value, double* parameters, double* action, int sync, int lmppx) {
  /*  auto& instance = instances.top();
    MPI_Comm ccomm = (lmppx != 0) ? MPI_COMM_SELF : instance.comm;
    int mpi_rank;
    MPI_Comm_rank(ccomm, &mpi_rank);
    Rvector ccc(instance.dimension, ccomm);
    auto ccrange = ccc.distribution().range(mpi_rank);
    auto ccn = ccrange.second - ccrange.first;
    ccc.allocate_buffer(Span<typename Rvector::value_type>(&parameters[ccrange.first], ccn));
    Rvector ggg(instance.dimension, ccomm);
    auto ggrange = ggg.distribution().range(mpi_rank);
    auto ggn = ggrange.second - ggrange.first;
    ggg.allocate_buffer(Span<typename Rvector::value_type>(&action[ggrange.first], ggn));
    size_t working_set_size =
        static_cast<Optimize<Rvector, Qvector>*>(instance.solver.get())->addValue(ccc, value, ggg) ? 1 : 0;
    if (sync) { // throw an error if communicator was not passed?
      gather_all(ccc.distribution(), ccomm, &parameters[0]);
      gather_all(ggg.distribution(), ccomm, &action[0]);
    }
    return working_set_size;*/
  return 0;
}

void apply_on_p_c(const std::vector<vectorP>& pvectors, const CVecRef<Pvector>& pspace, const VecRef<Rvector>& action) {
  auto& instance = instances.top();
  MPI_Comm ccomm =
      instance.comm; // TODO: not checking for MPI_COMM_SELF! Should lmppx be passed once and kept in instance?
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

extern "C" size_t IterativeSolverAddVector(double* parameters, double* action, int sync, int lmppx) {
  std::vector<Rvector> cc, gg;
  auto& instance = instances.top();
  if (instance.prof != nullptr)
    instance.prof->start("AddVector");
  cc.reserve(instance.solver->n_roots()); // TODO: should that be size of working set instead?
  gg.reserve(instance.solver->n_roots());
  MPI_Comm ccomm = (lmppx != 0) ? MPI_COMM_SELF : instance.comm;
  int mpi_rank;
  MPI_Comm_rank(ccomm, &mpi_rank);
  for (size_t root = 0; root < instance.solver->n_roots(); root++) {
    cc.emplace_back(instance.dimension, ccomm);
    auto ccrange = cc.back().distribution().range(mpi_rank);
    auto ccn = ccrange.second - ccrange.first;
    cc.back().allocate_buffer(
        Span<typename Rvector::value_type>(&parameters[root * instance.dimension + ccrange.first], ccn));
    gg.emplace_back(instance.dimension, ccomm);
    auto ggrange = gg.back().distribution().range(mpi_rank);
    auto ggn = ggrange.second - ggrange.first;
    gg.back().allocate_buffer(
        Span<typename Rvector::value_type>(&action[root * instance.dimension + ggrange.first], ggn));
  }
  if (instance.prof != nullptr)
    instance.prof->start("AddVector:Update");
  size_t working_set_size = instance.solver->add_vector(cc, gg);
  if (instance.prof != nullptr)
    instance.prof->stop("AddVector:Update");

  if (instance.prof != nullptr)
    instance.prof->start("AddVector:Sync");
  for (size_t root = 0; root < instance.solver->n_roots(); root++) {
    if (sync) {
      gather_all(cc[root].distribution(), ccomm, &parameters[root * instance.dimension]);
      gather_all(gg[root].distribution(), ccomm, &action[root * instance.dimension]);
    }
  }
  if (instance.prof != nullptr)
    instance.prof->stop("AddVector:Sync");
  if (mpi_rank == 0)
    instance.solver->report();
  if (instance.prof != nullptr)
    instance.prof->stop("AddVector");
  return working_set_size;
}

extern "C" size_t IterativeSolverPspaceAddVector(double* parameters, double* action, int sync, int lmppx,
                                                 void (*func)(const double*, double*, const size_t, const size_t*)) {
  std::vector<Rvector> cc, gg;
  auto& instance = instances.top();
  instance.apply_on_p_fort = func;
  if (instance.prof != nullptr)
    instance.prof->start("AddVector");
  cc.reserve(instance.solver->n_roots()); // TODO: should that be size of working set instead?
  gg.reserve(instance.solver->n_roots());
  MPI_Comm ccomm = (lmppx != 0) ? MPI_COMM_SELF : instance.comm;
  int mpi_rank;
  MPI_Comm_rank(ccomm, &mpi_rank);
  for (size_t root = 0; root < instance.solver->n_roots(); root++) {
    cc.emplace_back(instance.dimension, ccomm);
    auto ccrange = cc.back().distribution().range(mpi_rank);
    auto ccn = ccrange.second - ccrange.first;
    cc.back().allocate_buffer(
        Span<typename Rvector::value_type>(&parameters[root * instance.dimension + ccrange.first], ccn));
    gg.emplace_back(instance.dimension, ccomm);
    auto ggrange = gg.back().distribution().range(mpi_rank);
    auto ggn = ggrange.second - ggrange.first;
    gg.back().allocate_buffer(
        Span<typename Rvector::value_type>(&action[root * instance.dimension + ggrange.first], ggn));
  }
  using vectorP = std::vector<double>;
  using molpro::linalg::itsolv::CVecRef;
  using molpro::linalg::itsolv::VecRef;
  std::function<void(const std::vector<vectorP>&, const CVecRef<Pvector>&, const VecRef<Rvector>&)> apply_on_p =
      apply_on_p_c;
  if (instance.prof != nullptr)
    instance.prof->start("AddVector:Update");
  size_t working_set_size = instance.solver->add_vector(cc, gg, apply_on_p_c);
  if (instance.prof != nullptr)
    instance.prof->stop("AddVector:Update");

  if (instance.prof != nullptr)
    instance.prof->start("AddVector:Sync");
  for (size_t root = 0; root < instance.solver->n_roots(); root++) {
    if (sync) {
      gather_all(cc[root].distribution(), ccomm, &parameters[root * instance.dimension]);
      gather_all(gg[root].distribution(), ccomm, &action[root * instance.dimension]);
    }
  }
  if (instance.prof != nullptr)
    instance.prof->stop("AddVector:Sync");
  if (mpi_rank == 0)
    instance.solver->report();
  if (instance.prof != nullptr)
    instance.prof->stop("AddVector");
  return working_set_size;
}

extern "C" void IterativeSolverSolution(int nroot, int* roots, double* parameters, double* action, int sync,
                                        int lmppx) {
  std::vector<Rvector> cc, gg;
  auto& instance = instances.top();
  if (instance.prof != nullptr)
    instance.prof->start("Solution");
  cc.reserve(instance.solver->n_roots());
  gg.reserve(instance.solver->n_roots());
  MPI_Comm ccomm = (lmppx != 0) ? MPI_COMM_SELF : instance.comm;
  int mpi_rank;
  MPI_Comm_rank(ccomm, &mpi_rank);
  for (size_t root = 0; root < instance.solver->n_roots(); root++) {
    cc.emplace_back(instance.dimension, ccomm);
    auto ccrange = cc.back().distribution().range(mpi_rank);
    auto ccn = ccrange.second - ccrange.first;
    cc.back().allocate_buffer(
        Span<typename Rvector::value_type>(&parameters[root * instance.dimension + ccrange.first], ccn));
    gg.emplace_back(instance.dimension, ccomm);
    auto ggrange = gg.back().distribution().range(mpi_rank);
    auto ggn = ggrange.second - ggrange.first;
    gg.back().allocate_buffer(
        Span<typename Rvector::value_type>(&action[root * instance.dimension + ggrange.first], ggn));
  }
  const std::vector<int> croots(roots, roots + nroot);
  if (instance.prof != nullptr)
    instance.prof->start("Solution:Call");
  instance.solver->solution(croots, cc, gg);
  if (instance.prof != nullptr)
    instance.prof->stop("Solution:Call");

  if (instance.prof != nullptr)
    instance.prof->start("Solution:Sync");
  for (size_t root = 0; root < instance.solver->n_roots(); root++) {
    if (sync) {
      gather_all(cc[root].distribution(), ccomm, &parameters[root * instance.dimension]);
      gather_all(gg[root].distribution(), ccomm, &action[root * instance.dimension]);
    }
  }
  if (instance.prof != nullptr)
    instance.prof->stop("Solution:Sync");
  if (mpi_rank == 0)
    instance.solver->report();
  if (instance.prof != nullptr)
    instance.prof->stop("Solution");
}

extern "C" int IterativeSolverEndIteration(double* solution, double* residual, double* error, int sync, int lmppx) {
  std::vector<Rvector> cc, gg;
  auto& instance = instances.top();
  if (instance.prof != nullptr)
    instance.prof->start("EndIter");
  cc.reserve(instance.solver->n_roots());
  gg.reserve(instance.solver->n_roots());
  MPI_Comm ccomm = (lmppx != 0) ? MPI_COMM_SELF : instance.comm;
  int mpi_rank;
  MPI_Comm_rank(ccomm, &mpi_rank);
  for (size_t root = 0; root < instance.solver->n_roots(); root++) {
    cc.emplace_back(instance.dimension, ccomm);
    auto ccrange = cc.back().distribution().range(mpi_rank);
    auto ccn = ccrange.second - ccrange.first;
    cc.back().allocate_buffer(
        Span<typename Rvector::value_type>(&solution[root * instance.dimension + ccrange.first], ccn));
    gg.emplace_back(instance.dimension, ccomm);
    auto ggrange = gg.back().distribution().range(mpi_rank);
    auto ggn = ggrange.second - ggrange.first;
    gg.back().allocate_buffer(
        Span<typename Rvector::value_type>(&residual[root * instance.dimension + ggrange.first], ggn));
  }
  if (instance.prof != nullptr)
    instance.prof->start("EndIter:Call");
  int result = instance.solver->end_iteration(cc, gg);
  if (instance.prof != nullptr)
    instance.prof->stop("EndIter:Call");
  if (instance.prof != nullptr)
    instance.prof->start("AddVector:Sync");
  for (size_t root = 0; root < instance.solver->n_roots(); root++) {
    if (sync) {
      gather_all(cc[root].distribution(), ccomm, &solution[root * instance.dimension]);
      gather_all(gg[root].distribution(), ccomm, &residual[root * instance.dimension]);
    }
  }
  if (instance.prof != nullptr)
    instance.prof->stop("AddVector:Sync");
  for (size_t root = 0; root < instance.solver->n_roots(); root++) {
    error[root] = instance.solver->errors()[root];
  }
  if (instance.prof != nullptr)
    instance.prof->stop("EndIter");
  return result;
}

extern "C" size_t IterativeSolverAddP(size_t nP, const size_t* offsets, const size_t* indices,
                                      const double* coefficients, const double* pp, double* parameters, double* action,
                                      int sync, int lmppx,
                                      void (*func)(const double*, double*, const size_t, const size_t*)) {
  std::vector<Rvector> cc, gg;
  auto& instance = instances.top();
  instance.apply_on_p_fort = func;
  if (instance.prof != nullptr)
    instance.prof->start("AddP");
  cc.reserve(instance.solver->n_roots());
  gg.reserve(instance.solver->n_roots());
  MPI_Comm ccomm = (lmppx != 0) ? MPI_COMM_SELF : instance.comm;
  int mpi_rank;
  MPI_Comm_rank(ccomm, &mpi_rank);
  for (size_t root = 0; root < instance.solver->n_roots(); root++) {
    cc.emplace_back(instance.dimension, ccomm);
    auto ccrange = cc.back().distribution().range(mpi_rank);
    auto ccn = ccrange.second - ccrange.first;
    cc.back().allocate_buffer(Span<Rvector::value_type>(&parameters[root * instance.dimension + ccrange.first], ccn));
    gg.emplace_back(instance.dimension, ccomm);
    auto ggrange = gg.back().distribution().range(mpi_rank);
    auto ggn = ggrange.second - ggrange.first;
    gg.back().allocate_buffer(Span<Rvector::value_type>(&action[root * instance.dimension + ggrange.first], ggn));
  }
  std::vector<Pvector> Pvectors;
  Pvectors.reserve(nP);
  for (size_t p = 0; p < nP; p++) {
    // std::map<size_t, Rvector::value_type> ppp;
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
  for (size_t root = 0; root < instance.solver->n_roots(); root++) { //! TODO: should be working_set()?
    if (sync) {
      gather_all(cc[root].distribution(), ccomm, &parameters[root * instance.dimension]);
      gather_all(gg[root].distribution(), ccomm, &action[root * instance.dimension]);
    }
  }
  if (instance.prof != nullptr)
    instance.prof->stop("AddP:Sync");
  if (instance.prof != nullptr)
    instance.prof->stop("AddP");
  return working_set_size;
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
  LinearEigensystem<Rvector, Qvector, Pvector>* solver_cast =
      dynamic_cast<LinearEigensystem<Rvector, Qvector, Pvector>*>(instance.solver.get());
  if (solver_cast) {
    for (const auto& e : solver_cast->working_set_eigenvalues())
      eigenvalues[k++] = e;
  }
}

extern "C" size_t IterativeSolverSuggestP(const double* solution, const double* residual, size_t maximumNumber,
                                          double threshold, size_t* indices, int lmppx) {
  std::vector<Rvector> cc, gg;
  auto& instance = instances.top();
  if (instance.prof != nullptr)
    instance.prof->start("EndIter");
  cc.reserve(instance.solver->n_roots());
  gg.reserve(instance.solver->n_roots());
  MPI_Comm ccomm = (lmppx != 0) ? MPI_COMM_SELF : instance.comm;
  int mpi_rank;
  MPI_Comm_rank(ccomm, &mpi_rank);
  for (size_t root = 0; root < instance.solver->n_roots(); root++) {
    cc.emplace_back(instance.dimension, ccomm);
    auto ccrange = cc.back().distribution().range(mpi_rank);
    auto ccn = ccrange.second - ccrange.first;
    cc.back().allocate_buffer(
        Span<Rvector::value_type>(&const_cast<double*>(solution)[root * instance.dimension + ccrange.first], ccn));
    gg.emplace_back(instance.dimension, ccomm);
    auto ggrange = gg.back().distribution().range(mpi_rank);
    auto ggn = ggrange.second - ggrange.first;
    gg.back().allocate_buffer(
        Span<Rvector::value_type>(&const_cast<double*>(residual)[root * instance.dimension + ggrange.first], ggn));
  }
  auto result = instance.solver->suggest_p(molpro::linalg::itsolv::cwrap(cc), molpro::linalg::itsolv::cwrap(gg),
                                           maximumNumber, threshold);
  for (size_t i = 0; i < result.size(); i++) {
    indices[i] = result[i];
  }
  return result.size();
}

extern "C" void IterativeSolverPrintStatistics() { molpro::cout << instances.top().solver->statistics() << std::endl; }