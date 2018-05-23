#include "IterativeSolver.h"
#include "PagedVector.h"

// C interface to IterativeSolver
namespace LinearAlgebra {
using v = PagedVector<double>;

static std::unique_ptr<LinearEigensystem<double> > instance;

extern "C" void IterativeSolverLinearEigensystemInitialize(size_t n, size_t nroot, double thresh, int maxIterations, int verbosity) {
#ifdef USE_MPI
 int flag;
 MPI_Initialized(&flag);
 if (!flag) MPI_Init(0,nullptr);
#endif
 instance.reset(new LinearEigensystem<double>());
 instance->m_dimension=n;
 instance->m_roots=nroot;
 instance->m_thresh=thresh;
 instance->m_maxIterations=maxIterations;
 instance->m_verbosity=verbosity;
 std::cout << "roots: "<<instance->m_roots<<std::endl;
}
extern "C" void IterativeSolverLinearEigensystemAddVector(double* parameters, double* action, double* eigenvalue) {
 vectorSet<double> cc,gg;
 for (int root=0; root < instance->m_roots; root++) {
  cc.push_back(std::shared_ptr<v>(new v(instance->m_dimension)));
  cc.back()->put(&parameters[root*instance->m_dimension],instance->m_dimension,0);
  gg.push_back(std::shared_ptr<v>(new v(instance->m_dimension)));
  gg.back()->put(&action[root*instance->m_dimension],instance->m_dimension,0);
 }
 instance->addVector(cc,gg);
 for (int root=0; root < instance->m_roots; root++) {
  cc[root]->get(&parameters[root*instance->m_dimension],instance->m_dimension,0);
  gg[root]->get(&action[root*instance->m_dimension],instance->m_dimension,0);
  eigenvalue[root] = instance->eigenvalues()[root];
 }
}

extern "C" int IterativeSolverLinearEigensystemEndIteration(double* solution, double* residual, double* error) {
 vectorSet<double> cc,gg;
 for (int root=0; root < instance->m_roots; root++) {
  cc.push_back(std::shared_ptr<v>(new v(instance->m_dimension)));
  cc.back()->put(&solution[root*instance->m_dimension],instance->m_dimension,0);
  gg.push_back(std::shared_ptr<v>(new v(instance->m_dimension)));
  gg.back()->put(&residual[root*instance->m_dimension],instance->m_dimension,0);
 }
 bool result = instance->endIteration(cc,gg);
 for (int root=0; root < instance->m_roots; root++) {
  cc[root]->get(&solution[root*instance->m_dimension],instance->m_dimension,0);
  error[root] = instance->errors()[root];
 }
 return result;
}

extern "C" void IterativeSolverAddP(size_t* indices, double* coefficients, double* pp) {

}

}


