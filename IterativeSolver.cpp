#include "IterativeSolver.h"
#include "PagedVector.h"

// C interface to IterativeSolver
namespace LinearAlgebra {
 using v = PagedVector<double>;

 static std::unique_ptr<IterativeSolver<double> > instance;

 extern "C" void
 IterativeSolverLinearEigensystemInitialize(size_t n, size_t nroot, double thresh, int maxIterations, int verbosity) {
#ifdef USE_MPI
  int flag;
  MPI_Initialized(&flag);
  if (!flag) MPI_Init(0,nullptr);
#endif
  instance.reset(new LinearEigensystem<double>());
  instance->m_dimension = n;
  instance->m_roots = nroot;
  instance->m_thresh = thresh;
  instance->m_maxIterations = maxIterations;
  instance->m_verbosity = verbosity;
 }

 extern "C" void
 IterativeSolverDIISInitialize(size_t n, size_t nroot, double thresh, int maxIterations, int verbosity) {
#ifdef USE_MPI
  int flag;
  MPI_Initialized(&flag);
  if (!flag) MPI_Init(0,nullptr);
#endif
  instance.reset(new DIIS<double>());
  instance->m_dimension = n;
  instance->m_thresh = thresh;
  instance->m_maxIterations = maxIterations;
  instance->m_verbosity = verbosity;
 }

 extern "C" void IterativeSolverFinalize() {
  instance.release();
 }

 extern "C" void IterativeSolverAddVector(double *parameters, double *action, double *parametersP) {
  vectorSet<double> cc, gg;
  std::vector<std::vector<double> > ccp;
  for (int root = 0; root < instance->m_roots; root++) {
   cc.push_back(std::shared_ptr<v>(new v(instance->m_dimension)));
   cc.back()->put(&parameters[root * instance->m_dimension], instance->m_dimension, 0);
   gg.push_back(std::shared_ptr<v>(new v(instance->m_dimension)));
   gg.back()->put(&action[root * instance->m_dimension], instance->m_dimension, 0);
   cc.m_active[root]=gg.m_active[root]=instance->errors().size()<=root ||  instance->errors()[root]>=instance->m_thresh;
  }
  instance->addVector(cc, gg, ccp);
  for (int root = 0; root < instance->m_roots; root++) {
   cc[root]->get(&parameters[root * instance->m_dimension], instance->m_dimension, 0);
   gg[root]->get(&action[root * instance->m_dimension], instance->m_dimension, 0);
   for (size_t i = 0; i < ccp[0].size(); i++)
    parametersP[root * ccp[0].size() + i] = ccp[root][i];
  }
 }

 extern "C" int IterativeSolverEndIteration(double *solution, double *residual, double *error) {
  vectorSet<double> cc, gg;
  for (int root = 0; root < instance->m_roots; root++) {
   cc.push_back(std::shared_ptr<v>(new v(instance->m_dimension)));
   cc.back()->put(&solution[root * instance->m_dimension], instance->m_dimension, 0);
   gg.push_back(std::shared_ptr<v>(new v(instance->m_dimension)));
   gg.back()->put(&residual[root * instance->m_dimension], instance->m_dimension, 0);
  }
  bool result = instance->endIteration(cc, gg);
  for (int root = 0; root < instance->m_roots; root++) {
   cc[root]->get(&solution[root * instance->m_dimension], instance->m_dimension, 0);
   error[root] = instance->errors()[root];
  }
  return result;
 }

 extern "C" void IterativeSolverAddP(const size_t nP, const size_t *offsets, const size_t *indices,
                                                      const double *coefficients, const double *pp,
                                                      double *parameters, double *action, double *parametersP) {
  vectorSet<double> cc, gg;
  std::vector<std::vector<double> > ccp;
  for (int root = 0; root < instance->m_roots; root++) {
   cc.push_back(std::shared_ptr<v>(new v(instance->m_dimension)));
   gg.push_back(std::shared_ptr<v>(new v(instance->m_dimension)));
  }
  std::vector<std::map<size_t, double> > Pvectors;
  for (size_t p = 0; p < nP; p++) {
   std::map<size_t, double> pp;
   for (size_t k = offsets[p]; k < offsets[p + 1]; k++)
    pp.insert(std::pair<size_t, double>(indices[k], coefficients[k]));
   Pvectors.emplace_back(pp);
  }

  instance->addP(Pvectors, pp, cc, gg, ccp);
  for (int root = 0; root < instance->m_roots; root++) {
   cc[root]->get(&parameters[root * instance->m_dimension], instance->m_dimension, 0);
   gg[root]->get(&action[root * instance->m_dimension], instance->m_dimension, 0);
   for (size_t i = 0; i < ccp[0].size(); i++)
    parametersP[root * ccp[0].size() + i] = ccp[root][i];
  }
 }

 extern "C" void IterativeSolverEigenvalues(double* eigenvalues) {
  size_t k=0;
  for (const auto& e : instance->eigenvalues()) eigenvalues[k++] = e;
 }
}

