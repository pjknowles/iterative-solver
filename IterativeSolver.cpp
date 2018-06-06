#include "IterativeSolver.h"
#include "PagedVector.h"

// C interface to IterativeSolver
namespace LinearAlgebra {
 using v = PagedVector<double>;

 static std::unique_ptr<IterativeSolver<double> > instance;

 extern "C" void
 IterativeSolverLinearEigensystemInitialize(size_t n, size_t nroot, double thresh, unsigned int maxIterations, int verbosity) {
#ifdef HAVE_MPI_H
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
 IterativeSolverDIISInitialize(size_t n, double thresh, unsigned int maxIterations, int verbosity) {
#ifdef HAVE_MPI_H
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
//  constexpr int vFlags = LINEARALGEBRA_CLONE_ADVISE_DISTRIBUTED | LINEARALGEBRA_CLONE_ADVISE_OFFLINE; // has to wait till implementation below catches up
  constexpr int vFlags =  LINEARALGEBRA_CLONE_ADVISE_OFFLINE;
  vectorSet<double> cc, gg;
  std::vector<std::vector<double> > ccp;
  for (size_t root = 0; root < instance->m_roots; root++) {
   cc.push_back(std::shared_ptr<v>(new v(instance->m_dimension,vFlags)));
   cc.back()->put(&parameters[root * instance->m_dimension], instance->m_dimension, 0);
   gg.push_back(std::shared_ptr<v>(new v(instance->m_dimension,vFlags)));
   gg.back()->put(&action[root * instance->m_dimension], instance->m_dimension, 0);
   cc.m_active[root]=gg.m_active[root]=instance->errors().size()<=root ||  instance->errors()[root]>=instance->m_thresh;
  }
  instance->addVector(cc, gg, ccp);
  for (size_t root = 0; root < instance->m_roots; root++) {
   cc[root]->get(&parameters[root * instance->m_dimension], instance->m_dimension, 0);
   gg[root]->get(&action[root * instance->m_dimension], instance->m_dimension, 0);
#ifdef HAVE_MPI_H
   if (cc[root]->m_mpi_rank) throw std::logic_error("incomplete implementation");// FIXME not sure
#endif
   for (size_t i = 0; i < ccp[0].size(); i++)
    parametersP[root * ccp[0].size() + i] = ccp[root][i];
  }
 }

 extern "C" int IterativeSolverEndIteration(double *solution, double *residual, double *error) {
  vectorSet<double> cc, gg;
  for (size_t root = 0; root < instance->m_roots; root++) {
   cc.push_back(std::shared_ptr<v>(new v(instance->m_dimension)));
   cc.back()->put(&solution[root * instance->m_dimension], instance->m_dimension, 0);
   gg.push_back(std::shared_ptr<v>(new v(instance->m_dimension)));
   gg.back()->put(&residual[root * instance->m_dimension], instance->m_dimension, 0);
  }
  bool result = instance->endIteration(cc, gg);
  for (size_t root = 0; root < instance->m_roots; root++) {
   cc[root]->get(&solution[root * instance->m_dimension], instance->m_dimension, 0);
   error[root] = instance->errors()[root];
  }
#ifdef USE_MPI
  //   if (cc[root]->m_mpi_rank) throw std::logic_error("incomplete implementation");
#endif USE_MPI
  return result;
 }

 extern "C" void IterativeSolverAddP(const size_t nP, const size_t *offsets, const size_t *indices,
                                                      const double *coefficients, const double *pp,
                                                      double *parameters, double *action, double *parametersP) {
  vectorSet<double> cc, gg;
  std::vector<std::vector<double> > ccp;
  for (size_t root = 0; root < instance->m_roots; root++) {
   cc.push_back(std::shared_ptr<v>(new v(instance->m_dimension)));
   gg.push_back(std::shared_ptr<v>(new v(instance->m_dimension)));
  }
  std::vector<std::map<size_t, double> > Pvectors;
  for (size_t p = 0; p < nP; p++) {
   std::map<size_t, double> ppp;
   for (size_t k = offsets[p]; k < offsets[p + 1]; k++)
    ppp.insert(std::pair<size_t, double>(indices[k], coefficients[k]));
   Pvectors.emplace_back(ppp);
  }

  instance->addP(Pvectors, pp, cc, gg, ccp);
  for (size_t root = 0; root < instance->m_roots; root++) {
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
