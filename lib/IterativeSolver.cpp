#include "IterativeSolver.h"
//#include "PagedVector.h"
#include "OutOfCoreArray.h"
#include <memory>
#include <string>
#include <stack>

// C interface to IterativeSolver
//using v = LinearAlgebra::PagedVector<double>;
using v = LinearAlgebra::OutOfCoreArray<double>;

static std::stack<std::unique_ptr<IterativeSolver::Base<v> > > instances;

extern "C" void
IterativeSolverLinearEigensystemInitialize(size_t n,
                                           size_t nroot,
                                           double thresh,
                                           unsigned int maxIterations,
                                           int verbosity,
                                           int orthogonalize,
                                           int lmppx) {
#ifdef HAVE_MPI_H
  int flag;
  MPI_Initialized(&flag);
  MPI_Comm pcomm;
  if (!flag) {
     MPI_Init(0, nullptr);
     pcomm = MPI_COMM_WORLD;
  } else if (lmppx != 0) {
     pcomm = MPI_COMM_SELF;
  } else {
     pcomm = MPI_COMM_COMPUTE; //defined in OOCA header
  }
#endif
  instances.push(std::make_unique<IterativeSolver::LinearEigensystem<v> >(IterativeSolver::LinearEigensystem<v>()));
  auto& instance = instances.top();
  instance->m_dimension = n;
  instance->m_roots = nroot;
  instance->m_thresh = thresh;
  instance->m_maxIterations = maxIterations;
  instance->m_verbosity = verbosity;
  instance->m_orthogonalize = orthogonalize;
}

extern "C" void
IterativeSolverLinearEquationsInitialize(size_t n,
                                         size_t nroot,
                                         const double* rhs,
                                         double aughes,
                                         double thresh,
                                         unsigned int maxIterations,
                                         int verbosity,
                                         int orthogonalize) {
#ifdef HAVE_MPI_H
  int flag;
  MPI_Initialized(&flag);
  if (!flag) MPI_Init(0, nullptr);
#endif
  std::vector<v> rr;
  rr.reserve(nroot);
  for (size_t root = 0; root < nroot; root++) {
      rr.emplace_back(&const_cast<double*>(rhs)[root * n],n);
    //rr.push_back(v(const_cast<double*>(&rhs[root * n]),
    //               n)); // in principle the const_cast is dangerous, but we trust LinearEquations to behave
  }
    instances.push(std::make_unique<IterativeSolver::LinearEquations<v> >(rr,aughes));
  //instances.push(std::make_unique<IterativeSolver::LinearEquations<v> >(IterativeSolver::LinearEquations<v>(rr,
  //                                                                                                          aughes)));
  auto& instance = instances.top();
  instance->m_dimension = n;
  instance->m_roots = nroot;
  instance->m_thresh = thresh;
  instance->m_maxIterations = maxIterations;
  instance->m_verbosity = verbosity;
  instance->m_orthogonalize = orthogonalize;
}

extern "C" void
IterativeSolverDIISInitialize(size_t n, double thresh, unsigned int maxIterations, int verbosity) {
#ifdef HAVE_MPI_H
  int flag;
  MPI_Initialized(&flag);
  if (!flag) MPI_Init(0, nullptr);
#endif
  instances.push(std::make_unique<IterativeSolver::DIIS<v> >(IterativeSolver::DIIS<v>()));
  auto& instance = instances.top();
  instance->m_dimension = n;
  instance->m_thresh = thresh;
  instance->m_maxIterations = maxIterations;
  instance->m_verbosity = verbosity;
}
extern "C" void
IterativeSolverOptimizeInitialize(size_t n,
                                  double thresh,
                                  unsigned int maxIterations,
                                  int verbosity,
                                  char* algorithm,
                                  int minimize) {
#ifdef HAVE_MPI_H
  int flag;
  MPI_Initialized(&flag);
  if (!flag) MPI_Init(0, nullptr);
#endif
  if (*algorithm)
    instances.push(std::make_unique<IterativeSolver::Optimize<v> >(IterativeSolver::Optimize<v>(algorithm,
                                                                                                minimize != 0)));
  else
    instances.push(std::make_unique<IterativeSolver::Optimize<v> >(IterativeSolver::Optimize<v>()));
  auto& instance = instances.top();
  instance->m_dimension = n;
  instance->m_roots = 1;
  instance->m_thresh = thresh;
  instance->m_maxIterations = maxIterations;
  instance->m_verbosity = verbosity;
}

extern "C" void IterativeSolverFinalize() {
  instances.pop();
}

extern "C" int IterativeSolverAddValue(double* parameters, double value, double* action, int sync, int lmppx) {
  auto& instance = instances.top();
#ifdef HAVE_MPI_H
  MPI_Comm ccomm;
  if (lmppx != 0) { //OK?
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
  auto result = static_cast<IterativeSolver::Optimize<v>*>(instance.get())->addValue(ccc, value, ggg) ? 1 : 0;
#ifdef HAVE_MPI_H
  if (sync) { //throw an error if communicator was not passed?
    if (!ccc.synchronised()) ccc.sync();
    if (!ggg.synchronised()) ggg.sync();
  }
#endif
  return result;
}

extern "C" int IterativeSolverAddVector(double* parameters, double* action, double* parametersP, int sync, int lmppx) {
  std::vector<v> cc, gg;
  auto& instance = instances.top();
  cc.reserve(instance->m_roots); // very important for avoiding copying of memory-mapped vectors in emplace_back below
  gg.reserve(instance->m_roots);
  std::vector<std::vector<typename v::value_type> > ccp(instance->m_roots);
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
  bool update = instance->addVector(cc, gg, ccp);

  for (size_t root = 0; root < instance->m_roots; root++) {
#ifdef HAVE_MPI_H
    if (sync) {
      if (!cc[root].synchronised()) cc[root].sync();
      if (!gg[root].synchronised()) gg[root].sync();
    }
#endif
    for (size_t i = 0; i < ccp[0].size(); i++)
      parametersP[root * ccp[0].size() + i] = ccp[root][i];
  }
  return update ? 1 : 0;
}

extern "C" int IterativeSolverEndIteration(double* solution, double* residual, double* error, int lmppx) {
  std::vector<v> cc, gg;
  auto& instance = instances.top();
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
  bool result = instance->endIteration(cc, gg);
  for (size_t root = 0; root < instance->m_roots; root++) {
#ifdef HAVE_MPI_H
    if (!cc[root].synchronised()) cc[root].sync();
    if (!gg[root].synchronised()) gg[root].sync();
#endif
    error[root] = instance->errors()[root];
  }
  return result;
}

extern "C" void IterativeSolverAddP(size_t nP, const size_t* offsets, const size_t* indices,
                                    const double* coefficients, const double* pp,
                                    double* parameters, double* action, double* parametersP,
                                    int lmppx) {
  std::vector<v> cc, gg;
  auto& instance = instances.top();
  std::vector<std::vector<v::value_type> > ccp(instance->m_roots);
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
    cc.emplace_back(&parameters[root * instance->m_dimension],instance->m_dimension, ccomm);
    gg.emplace_back(&action[root * instance->m_dimension],instance->m_dimension, ccomm);
  }
#else
  for (size_t root = 0; root < instance->m_roots; root++) {
    cc.emplace_back(&parameters[root * instance->m_dimension],instance->m_dimension);
    gg.emplace_back(&action[root * instance->m_dimension],instance->m_dimension);
  }
#endif
  std::vector<std::map<size_t, v::value_type> > Pvectors;
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
    if (!cc[root].synchronised()) cc[root].sync();
    if (!gg[root].synchronised()) gg[root].sync();
#endif
//    cc[root].get(&parameters[root * instance->m_dimension], instance->m_dimension, 0);
//    gg[root].get(&action[root * instance->m_dimension], instance->m_dimension, 0);
    for (size_t i = 0; i < ccp[0].size(); i++)
      parametersP[root * ccp[0].size() + i] = ccp[root][i];
  }
}

extern "C" void IterativeSolverOption(const char* key, const char* val) {
  auto& instance = instances.top();
  instance->m_options.insert(std::make_pair(std::string(key), std::string(val)));
}

extern "C" void IterativeSolverEigenvalues(double* eigenvalues) {
  auto& instance = instances.top();
  size_t k = 0;
  for (const auto& e : instance->eigenvalues()) eigenvalues[k++] = e;
}

extern "C" size_t IterativeSolverSuggestP(const double* solution,
                                          const double* residual,
                                          size_t maximumNumber,
                                          double threshold,
                                          size_t* indices,
                                          int lmppx) {
  std::vector<v> cc, gg;
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
    cc.emplace_back(&const_cast<double*>(solution)[root * instance->m_dimension],
                   instance->m_dimension, ccomm);
    gg.emplace_back(&const_cast<double*>(residual)[root * instance->m_dimension],
                   instance->m_dimension, ccomm);
  }
#else
  for (size_t root = 0; root < instance->m_roots; root++) {
    cc.emplace_back(&const_cast<double*>(solution)[root * instance->m_dimension],
                   instance->m_dimension);
    gg.emplace_back(&const_cast<double*>(residual)[root * instance->m_dimension],
                   instance->m_dimension);
  }
#endif

  auto result = instance->suggestP(cc, gg, maximumNumber, threshold);
  for (size_t i = 0; i < result.size(); i++)
    indices[i] = result[i];
  return result.size();
}
