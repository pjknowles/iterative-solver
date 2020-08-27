#include "IterativeSolverC.h"
#include "IterativeSolver.h"
#include "OutOfCoreArray.h"
#include "molpro/ProfilerSingle.h"
#include <memory>
#include <stack>
#include <string>
#include <tuple>
#include <mpi.h>
#ifdef HAVE_PPIDD_H
#include <ppidd.h>
#endif

#include <molpro/linalg/array/DistrArrayMPI3.h>
#include <molpro/linalg/iterativesolver/ArrayHandlers.h>
#include <molpro/linalg/array/util/Distribution.h>
#include <molpro/linalg/array/Span.h>
#include <molpro/linalg/array/util/gather_all.h>

using molpro::linalg::DIIS;
using molpro::linalg::LinearEigensystem;
using molpro::linalg::LinearEquations;
using molpro::linalg::Optimize;
using molpro::linalg::iterativesolver::ArrayHandlers;
using molpro::linalg::array::Span;
using molpro::linalg::array::util::gather_all;

using Rvector = molpro::linalg::array::DistrArrayMPI3;
using Qvector = molpro::linalg::array::DistrArrayMPI3;
using Pvector = std::map<size_t, double>;

MPI_Comm commun;

// FIXME Only top solver is active at any one time. This should be documented somewhere.
static std::stack<std::unique_ptr<molpro::linalg::IterativeSolver<Rvector, Qvector, Pvector>>> instances;

extern "C" void IterativeSolverLinearEigensystemInitialize(size_t n, size_t nroot, size_t range_begin, size_t range_end,
                                                           double thresh, int verbosity, const char* fname,
                                                           int64_t fcomm, int lmppx) {
  std::shared_ptr<molpro::Profiler> profiler = nullptr;
  std::string pname(fname);
  int flag;
  MPI_Initialized(&flag);
  MPI_Comm pcomm;
  if (!flag) {
#ifdef HAVE_PPIDD_H
    PPIDD_Initialize(0, nullptr, PPIDD_IMPL_DEFAULT);
    pcomm = MPI_Comm_f2c(PPIDD_Worker_comm());
    commun = MPI_Comm_f2c(PPIDD_Worker_comm());
#else
    MPI_Init(0, nullptr);
    pcomm = MPI_COMM_WORLD;
    commun = MPI_COMM_WORLD;
#endif
  } else if (lmppx != 0) {
    pcomm = MPI_COMM_SELF;
    commun = MPI_COMM_SELF;
  } else {
    // TODO: Check this is safe. Will crash if handle is invalid.
    pcomm = MPI_Comm_f2c(fcomm);
    commun = MPI_Comm_f2c(fcomm);
  }
  // TODO: what if lmppx != 0 ?
  if (!pname.empty()) {
    profiler = molpro::ProfilerSingle::instance(pname, pcomm);
  }
  int mpi_rank;
  MPI_Comm_rank(commun, &mpi_rank);
  auto handlers = std::make_shared<ArrayHandlers<Rvector, Qvector, Pvector>>();
  /*instances.push(std::make_unique<LinearEigensystem<Rvector, Qvector, Pvector>>(
      LinearEigensystem<Rvector, Qvector, Pvector>(handlers, profiler)));*/
  instances.push(std::make_unique<LinearEigensystem<Rvector, Qvector, Pvector>>(handlers, profiler));
  auto& instance = instances.top();
  instance->m_dimension = n;
  instance->m_roots = nroot;
  instance->m_thresh = thresh;
  instance->m_verbosity = verbosity;
  std::vector<Rvector> x;
  std::vector<Rvector> g;
  for (size_t root = 0; root < instance->m_roots; root++) {
    x.emplace_back(n, commun);
    g.emplace_back(n, commun);
  }
  std::tie(range_begin, range_end) = x[0].distribution().range(mpi_rank);
}
  
extern "C" void IterativeSolverLinearEquationsInitialize(size_t n, size_t nroot, const double* rhs, double aughes,
                                                         double thresh, unsigned int maxIterations, int verbosity) {
/*#ifdef HAVE_MPI_H
  int flag;
  MPI_Initialized(&flag);
#ifdef HAVE_PPIDD_H
  if (!flag)
    PPIDD_Initialize(0, nullptr, PPIDD_IMPL_DEFAULT);
#else
  if (!flag)
    MPI_Init(0, nullptr);
#endif
#endif
  std::vector<Rvector> rr;
  rr.reserve(nroot);
  for (size_t root = 0; root < nroot; root++) {
    rr.emplace_back(n, MPI_COMM_COMPUTE);
    rr.back().allocate_buffer(Span<typename Rvector::double>(&const_cast<double*>(rhs)[root * n], n));
    // rr.push_back(v(const_cast<double*>(&rhs[root * n]),
    //               n)); // in principle the const_cast is dangerous, but we trust LinearEquations to behave
  }
  instances.push(std::make_unique<LinearEquations<Rvector, Qvector, Pvector>>(rr,
                                                              ArrayHandlers<Rvector, Qvector, Pvector>{}, aughes));
  // instances.push(std::make_unique<IterativeSolver::LinearEquations<v> >(IterativeSolver::LinearEquations<v>(rr,
  //                                                                                                          aughes)));
  auto& instance = instances.top();
  instance->m_dimension = n;
  instance->m_roots = nroot;
  instance->m_thresh = thresh;
  instance->m_maxIterations = maxIterations;
  instance->m_verbosity = verbosity;*/
}

extern "C" void IterativeSolverDIISInitialize(size_t n, double thresh, unsigned int maxIterations, int verbosity) {
/*#ifdef HAVE_MPI_H
  int flag;
  MPI_Initialized(&flag);
#ifdef HAVE_PPIDD_H
  if (!flag)
    PPIDD_Initialize(0, nullptr, PPIDD_IMPL_DEFAULT);
#else
  if (!flag)
    MPI_Init(0, nullptr);
#endif
#endif
  instances.emplace(std::make_unique<DIIS<v>>(make_handlers()));
  auto& instance = instances.top();
  instance->m_dimension = n;
  instance->m_thresh = thresh;
  instance->m_maxIterations = maxIterations;
  instance->m_verbosity = verbosity;*/
}
extern "C" void IterativeSolverOptimizeInitialize(size_t n, double thresh, unsigned int maxIterations, int verbosity,
                                                  char* algorithm, int minimize) {
/*#ifdef HAVE_MPI_H
  int flag;
  MPI_Initialized(&flag);
#ifdef HAVE_PPIDD_H
  if (!flag)
    PPIDD_Initialize(0, nullptr, PPIDD_IMPL_DEFAULT);
#else
  if (!flag)
    MPI_Init(0, nullptr);
#endif
#endif
  if (*algorithm)
    instances.emplace(std::make_unique<Optimize<v>>(make_handlers(), algorithm, minimize != 0));
  else
    instances.emplace(std::make_unique<Optimize<v>>(make_handlers()));
  auto& instance = instances.top();
  instance->m_dimension = n;
  instance->m_roots = 1;
  instance->m_thresh = thresh;
  instance->m_maxIterations = maxIterations;
  instance->m_verbosity = verbosity;*/
}

extern "C" void IterativeSolverFinalize() { instances.pop(); }

extern "C" int IterativeSolverAddValue(double* parameters, double value, double* action, int sync, int lmppx) {
/*  auto& instance = instances.top();
#ifdef HAVE_MPI_H
  MPI_Comm ccomm;
  if (lmppx != 0) { // OK?
    ccomm = MPI_COMM_SELF;
  } else {
    ccomm = MPI_COMM_COMPUTE;
  }
  v ccc(parameters, instance->m_dimension, ccomm);
  v ggg(action, instance->m_dimension, ccomm);
#else
  v ccc(parameters, instance->m_dimension);
  v ggg(action, instance->m_dimension);
#endif
  auto result = static_cast<molpro::linalg::Optimize<v>*>(instance.get())->addValue(ccc, value, ggg) ? 1 : 0;
#ifdef HAVE_MPI_H
  if (sync) { // throw an error if communicator was not passed?
    if (!ccc.synchronised())
      ccc.sync();
    if (!ggg.synchronised())
      ggg.sync();
  }
#endif
  return result;*/
  return 0;
}

extern "C" int IterativeSolverAddVector(double* parameters, double* action, double* parametersP, int sync, int lmppx) {
  std::vector<Rvector> cc, gg;
  auto& instance = instances.top();
  if (instance->m_profiler != nullptr)
    instance->m_profiler->start("AddVector");
  cc.reserve(instance->m_roots); // TODO: should that be size of working set instead?
  gg.reserve(instance->m_roots);
  std::vector<std::vector<typename Rvector::value_type>> ccp(instance->m_roots);
  MPI_Comm ccomm;
  if (lmppx != 0) {
    ccomm = MPI_COMM_SELF;
  } else {
    ccomm = commun;
  }
  int mpi_rank;
  MPI_Comm_rank(ccomm, &mpi_rank);
  for (size_t root = 0; root < instance->m_roots; root++) {
    cc.emplace_back(instance->m_dimension, ccomm);
    auto ccrange = cc.back().distribution().range(mpi_rank);
    auto ccn = ccrange.second - ccrange.first;
    cc.back().allocate_buffer(Span<typename Rvector::value_type>(&parameters[root * instance->m_dimension +
                                                                                                  ccrange.first], ccn));
    gg.emplace_back(instance->m_dimension, ccomm);
    auto ggrange = gg.back().distribution().range(mpi_rank);
    auto ggn = ggrange.second - ggrange.first;
    gg.back().allocate_buffer(Span<typename Rvector::value_type>(&action[root * instance->m_dimension +
                                                                                                  ggrange.first], ggn));
  }
  if (instance->m_profiler != nullptr)
    instance->m_profiler->start("AddVector:Update");
  size_t working_set_size = instance->addVector(cc, gg, ccp);
  if (instance->m_profiler != nullptr)
    instance->m_profiler->stop("AddVector:Update");

  if (instance->m_profiler != nullptr)
    instance->m_profiler->start("AddVector:Sync");
  for (size_t root = 0; root < instance->m_roots; root++) {
    if (sync) {
      gather_all(cc[root].distribution(), ccomm, &parameters[root * instance->m_dimension]);
      gather_all(gg[root].distribution(), ccomm, &action[root * instance->m_dimension]);
    }
    for (size_t i = 0; i < ccp[0].size(); i++)
      parametersP[root * ccp[0].size() + i] = ccp[root][i];
  }
  if (instance->m_profiler != nullptr)
    instance->m_profiler->stop("AddVector:Sync");
  if (mpi_rank == 0) instance->report();
  if (instance->m_profiler != nullptr)
    instance->m_profiler->stop("AddVector");
  return working_set_size;
}

extern "C" void IterativeSolverSolution(int nroot, int* roots, double* parameters, double* action, double* parametersP,
                                                                                                  int sync, int lmppx) {
  std::vector<Rvector> cc, gg;
  auto& instance = instances.top();
  if (instance->m_profiler != nullptr)
    instance->m_profiler->start("Solution");
  cc.reserve(instance->m_roots);
  gg.reserve(instance->m_roots);
  std::vector<std::vector<typename Rvector::value_type>> ccp(instance->m_roots);
  MPI_Comm ccomm;
  if (lmppx != 0) {
    ccomm = MPI_COMM_SELF;
  } else {
    ccomm = commun;
  }
  int mpi_rank;
  MPI_Comm_rank(ccomm, &mpi_rank);
  for (size_t root = 0; root < instance->m_roots; root++) {
    cc.emplace_back(instance->m_dimension, ccomm);
    auto ccrange = cc.back().distribution().range(mpi_rank);
    auto ccn = ccrange.second - ccrange.first;
    cc.back().allocate_buffer(Span<typename Rvector::value_type>(&parameters[root * instance->m_dimension +
                                                                             ccrange.first], ccn));
    gg.emplace_back(instance->m_dimension, ccomm);
    auto ggrange = gg.back().distribution().range(mpi_rank);
    auto ggn = ggrange.second - ggrange.first;
    gg.back().allocate_buffer(Span<typename Rvector::value_type>(&action[root * instance->m_dimension +
                                                                         ggrange.first], ggn));
  }
  const std::vector<int> croots(roots, roots+nroot);
  if (instance->m_profiler != nullptr)
    instance->m_profiler->start("Solution:Call");
  instance->solution(croots, cc, gg, ccp);
  if (instance->m_profiler != nullptr)
    instance->m_profiler->stop("Solution:Call");

  if (instance->m_profiler != nullptr)
    instance->m_profiler->start("Solution:Sync");
  for (size_t root = 0; root < instance->m_roots; root++) {
    if (sync) {
      gather_all(cc[root].distribution(), ccomm, &parameters[root * instance->m_dimension]);
      gather_all(gg[root].distribution(), ccomm, &action[root * instance->m_dimension]);
    }
    for (size_t i = 0; i < ccp[0].size(); i++)
      parametersP[root * ccp[0].size() + i] = ccp[root][i];
  }
  if (instance->m_profiler != nullptr)
    instance->m_profiler->stop("Solution:Sync");
  if (mpi_rank == 0) instance->report();
  if (instance->m_profiler != nullptr)
    instance->m_profiler->stop("Solution");
}

extern "C" int IterativeSolverEndIteration(double* solution, double* residual, double* error, int lmppx) {
/*  std::vector<v> cc, gg;
  auto& instance = instances.top();
  if (instance->m_profiler != nullptr)
    instance->m_profiler->start("EndIter");
  cc.reserve(instance->m_roots); // very important for avoiding copying of memory-mapped vectors in emplace_back below
  gg.reserve(instance->m_roots);
#ifdef HAVE_MPI_H
  MPI_Comm ccomm;
  if (lmppx != 0) {
    ccomm = MPI_COMM_SELF;
  } else {
    ccomm = MPI_COMM_COMPUTE;
  }
  for (size_t root = 0; root < instance->m_roots; root++) {
    cc.emplace_back(&solution[root * instance->m_dimension], instance->m_dimension, ccomm);
    gg.emplace_back(&residual[root * instance->m_dimension], instance->m_dimension, ccomm);
  }
#else
  for (size_t root = 0; root < instance->m_roots; root++) {
    cc.emplace_back(&solution[root * instance->m_dimension], instance->m_dimension);
    gg.emplace_back(&residual[root * instance->m_dimension], instance->m_dimension);
  }
#endif
  if (instance->m_profiler != nullptr)
    instance->m_profiler->start("EndIter:Body");
  bool result = instance->endIteration(cc, gg);
  if (instance->m_profiler != nullptr)
    instance->m_profiler->stop("EndIter:Body");
  if (instance->m_profiler != nullptr)
    instance->m_profiler->start("EndIter:Sync");
  for (size_t root = 0; root < instance->m_roots; root++) {
#ifdef HAVE_MPI_H
    if (!cc[root].synchronised())
      cc[root].sync();
    if (!gg[root].synchronised())
      gg[root].sync();
#endif
    error[root] = instance->errors()[root];
  }
  if (instance->m_profiler != nullptr)
    instance->m_profiler->stop("EndIter:Sync");
  if (instance->m_profiler != nullptr)
    instance->m_profiler->stop("EndIter");
  return result;*/
  return 0;
}

extern "C" void IterativeSolverAddP(size_t nP, const size_t* offsets, const size_t* indices, const double* coefficients,
                                    const double* pp, double* parameters, double* action, double* parametersP,
                                    int lmppx) {
/*  std::vector<v> cc, gg;
  auto& instance = instances.top();
  std::vector<std::vector<v::value_type>> ccp(instance->m_roots);
  cc.reserve(instance->m_roots); // very important for avoiding copying of memory-mapped vectors in emplace_back below
  gg.reserve(instance->m_roots);
#ifdef HAVE_MPI_H
  MPI_Comm ccomm;
  if (lmppx != 0) {
    ccomm = MPI_COMM_SELF;
  } else {
    ccomm = MPI_COMM_COMPUTE;
  }
  for (size_t root = 0; root < instance->m_roots; root++) {
    cc.emplace_back(&parameters[root * instance->m_dimension], instance->m_dimension, ccomm);
    gg.emplace_back(&action[root * instance->m_dimension], instance->m_dimension, ccomm);
  }
#else
  for (size_t root = 0; root < instance->m_roots; root++) {
    cc.emplace_back(&parameters[root * instance->m_dimension], instance->m_dimension);
    gg.emplace_back(&action[root * instance->m_dimension], instance->m_dimension);
  }
#endif
  std::vector<std::map<size_t, v::value_type>> Pvectors;
  Pvectors.reserve(nP);
  for (size_t p = 0; p < nP; p++) {
    std::map<size_t, v::value_type> ppp;
    //    for (size_t k = offsets[p]; k < offsets[p + 1]; k++)
    //    std::cout << "indices["<<k<<"]="<<indices[k]<<": "<<coefficients[k]<<std::endl;
    for (size_t k = offsets[p]; k < offsets[p + 1]; k++)
      ppp.insert(std::pair<size_t, v::value_type>(indices[k], coefficients[k]));
    Pvectors.emplace_back(ppp);
  }

  instance->addP(Pvectors, pp, cc, gg, ccp);
  for (size_t root = 0; root < instance->m_roots; root++) {
#ifdef HAVE_MPI_H
    if (!cc[root].synchronised())
      cc[root].sync();
    if (!gg[root].synchronised())
      gg[root].sync();
#endif
    //    cc[root].get(&parameters[root * instance->m_dimension], instance->m_dimension, 0);
    //    gg[root].get(&action[root * instance->m_dimension], instance->m_dimension, 0);
    for (size_t i = 0; i < ccp[0].size(); i++)
      parametersP[root * ccp[0].size() + i] = ccp[root][i];
  }*/
}

extern "C" void IterativeSolverEigenvalues(double* eigenvalues) {
  auto& instance = instances.top();
  size_t k = 0;
  for (const auto& e : instance->eigenvalues())
    eigenvalues[k++] = e;
}

extern "C" void IterativeSolverWorkingSetEigenvalues(double* eigenvalues) {
  auto& instance = instances.top();
  size_t k = 0;
  for (const auto& e : instance->working_set_eigenvalues())
    eigenvalues[k++] = e;
}

extern "C" size_t IterativeSolverSuggestP(const double* solution, const double* residual, size_t maximumNumber,
                                          double threshold, size_t* indices, int lmppx) {
/*  std::vector<v> cc, gg;
  auto& instance = instances.top();
  cc.reserve(instance->m_roots);
  gg.reserve(instance->m_roots);
#ifdef HAVE_MPI_H
  MPI_Comm ccomm;
  if (lmppx != 0) {
    ccomm = MPI_COMM_SELF;
  } else {
    ccomm = MPI_COMM_COMPUTE;
  }
  for (size_t root = 0; root < instance->m_roots; root++) {
    cc.emplace_back(&const_cast<double*>(solution)[root * instance->m_dimension], instance->m_dimension, ccomm);
    gg.emplace_back(&const_cast<double*>(residual)[root * instance->m_dimension], instance->m_dimension, ccomm);
  }
#else
  for (size_t root = 0; root < instance->m_roots; root++) {
    cc.emplace_back(&const_cast<double*>(solution)[root * instance->m_dimension], instance->m_dimension);
    gg.emplace_back(&const_cast<double*>(residual)[root * instance->m_dimension], instance->m_dimension);
  }
#endif

  auto result = instance->suggestP(cc, gg, maximumNumber, threshold);
  for (size_t i = 0; i < result.size(); i++)
    indices[i] = result[i];
  return result.size();*/
  return 0;
}
