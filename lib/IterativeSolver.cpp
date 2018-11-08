#include "IterativeSolver.h"
#include "PagedVector.h"

//namespace LinearAlgebra {
//template
//class LinearEigensystem<double>;
//
//template
//class LinearEquations<double>;
//
//template
//class DIIS<double>;
//}

// C interface to IterativeSolver
namespace LinearAlgebra {
using v = PagedVector<double>;

static std::unique_ptr<IterativeSolver<v> > instance;

extern "C" void
IterativeSolverLinearEigensystemInitialize(size_t n,
                                           size_t nroot,
                                           double thresh,
                                           unsigned int maxIterations,
                                           int verbosity,
                                           int orthogonalize) {
#ifdef HAVE_MPI_H
  int flag;
  MPI_Initialized(&flag);
  if (!flag) MPI_Init(0, nullptr);
#endif
  instance.reset(new LinearEigensystem<v>());
  instance->m_dimension = n;
  instance->m_roots = nroot;
  instance->m_thresh = thresh;
  instance->m_maxIterations = maxIterations;
  instance->m_verbosity = verbosity;
  instance->m_orthogonalize = orthogonalize;
  std::cout << "orthogonalize: " << orthogonalize << std::endl;
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
  for (size_t root = 0; root < nroot; root++) {
    rr.push_back(v(const_cast<double*>(&rhs[root * n]),
                   n)); // in principle the const_cast is dangerous, but we trust LinearEquations to behanve
  }
  instance.reset(new LinearEquations<v>(rr, aughes));
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
  instance.reset(new DIIS<v>());
  instance->m_dimension = n;
  instance->m_thresh = thresh;
  instance->m_maxIterations = maxIterations;
  instance->m_verbosity = verbosity;
}

extern "C" void IterativeSolverFinalize() {
  instance.release();
}

extern "C" void IterativeSolverAddVector(double* parameters, double* action, const int* active, double* parametersP) {
//  constexpr int vFlags = LINEARALGEBRA_CLONE_ADVISE_DISTRIBUTED | LINEARALGEBRA_CLONE_ADVISE_OFFLINE; // has to wait till implementation below catches up
//  constexpr int vFlags = LINEARALGEBRA_CLONE_ADVISE_OFFLINE; // FIXME is this needed?
  std::vector<v> cc, gg;
  cc.reserve(instance->m_roots); // very important for avoiding copying of memory-mapped vectors in emplace_back below
  gg.reserve(instance->m_roots);
  std::vector<std::vector<typename v::value_type> > ccp;
  std::vector<bool> activev;
  for (size_t root = 0; root < instance->m_roots; root++) {
    cc.emplace_back(&parameters[root * instance->m_dimension], instance->m_dimension);
    gg.emplace_back(&action[root * instance->m_dimension], instance->m_dimension);
//    activev.push_back( instance->errors().size() <= root || instance->errors()[root] >= instance->m_thresh);
    activev.push_back(active[root]!=0);
  }
  instance->addVector(cc, gg, activev, ccp);

  for (size_t root = 0; root < instance->m_roots; root++) {
    for (size_t i = 0; i < ccp[0].size(); i++)
      parametersP[root * ccp[0].size() + i] = ccp[root][i];
  }
}

extern "C" int IterativeSolverEndIteration(double* solution, double* residual, double* error, int* active) {
  std::vector<v> cc, gg;
  cc.reserve(instance->m_roots); // very important for avoiding copying of memory-mapped vectors in emplace_back below
  gg.reserve(instance->m_roots);
  std::vector<bool> activev;
  for (size_t root = 0; root < instance->m_roots; root++) {
    activev.push_back(1);
    cc.emplace_back(&solution[root * instance->m_dimension], instance->m_dimension);
    gg.emplace_back(&residual[root * instance->m_dimension], instance->m_dimension);
  }
  bool result = instance->endIteration(cc, gg, activev);
  for (size_t root = 0; root < instance->m_roots; root++) {
//    cc[root].get(&solution[root * instance->m_dimension], instance->m_dimension, 0);
    error[root] = instance->errors()[root];
    active[root] = activev[root] ? 1 : 0;
  }
#ifdef HAVE_MPI_H
  //   if (cc[root]->m_mpi_rank) throw std::logic_error("incomplete implementation");
#endif
  return result;
}

extern "C" void IterativeSolverAddP(size_t nP, const size_t* offsets, const size_t* indices,
                                    const double* coefficients, const double* pp,
                                    double* parameters, double* action, double* parametersP) {
  std::vector<v> cc, gg;
  std::vector<std::vector<v::value_type> > ccp;
  for (size_t root = 0; root < instance->m_roots; root++) {
    cc.push_back(v(instance->m_dimension));
    gg.push_back(v(instance->m_dimension));
  }
  std::vector<std::map<size_t, v::value_type> > Pvectors;
  Pvectors.reserve(nP);
  for (size_t p = 0; p < nP; p++) {
    std::map<size_t, v::value_type> ppp;
    for (size_t k = offsets[p]; k < offsets[p + 1]; k++)
//    std::cout << "indices["<<k<<"]="<<indices[k]<<": "<<coefficients[k]<<std::endl;
      for (size_t k = offsets[p]; k < offsets[p + 1]; k++)
        ppp.insert(std::pair<size_t, v::value_type>(indices[k], coefficients[k]));
    Pvectors.emplace_back(ppp);
  }

  instance->addP(Pvectors, pp, cc, gg, ccp);
  for (size_t root = 0; root < instance->m_roots; root++) {
    cc[root].get(&parameters[root * instance->m_dimension], instance->m_dimension, 0);
    gg[root].get(&action[root * instance->m_dimension], instance->m_dimension, 0);
    for (size_t i = 0; i < ccp[0].size(); i++)
      parametersP[root * ccp[0].size() + i] = ccp[root][i];
  }
}

extern "C" void IterativeSolverOption(const char* key, const char* val) {
  instance->m_options.insert(std::make_pair(std::string(key), std::string(val)));
}

extern "C" void IterativeSolverEigenvalues(double* eigenvalues) {
  size_t k = 0;
  for (const auto& e : instance->eigenvalues()) eigenvalues[k++] = e;
}

extern "C" size_t IterativeSolverSuggestP(const double* solution,
                                          const double* residual,
                                          const int* activec,
                                          size_t maximumNumber,
                                          double threshold,
                                          size_t* indices) {
  std::vector<v> cc, gg;
  cc.reserve(instance->m_roots);
  gg.reserve(instance->m_roots);
  std::vector<bool> active;
  for (size_t root = 0; root < instance->m_roots; root++) {
    active.push_back(activec[root]!=0);
    cc.push_back(v(&const_cast<double*>(solution)[root * instance->m_dimension],
                   instance->m_dimension));
    gg.push_back(v(&const_cast<double*>(residual)[root * instance->m_dimension],
                   instance->m_dimension));
    active[root] = (instance->errors().size() <= root || instance->errors()[root] >= instance->m_thresh);
  }

  auto result = instance->suggestP(cc, gg, active, maximumNumber, threshold);
  for (size_t i = 0; i < result.size(); i++)
    indices[i] = result[i];
  return result.size();
}
}
