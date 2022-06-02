#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_ITERATIVESOLVERC_H_
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_ITERATIVESOLVERC_H_
#include <cstdint>
#include <stddef.h>

extern "C" void IterativeSolverLinearEigensystemInitialize(size_t nQ, size_t nroot, size_t* range_begin,
                                                           size_t* range_end, double thresh, double thresh_value,
                                                           int hermitian, int verbosity, const char* fname,
                                                           int64_t fcomm, const char* algorithm, const char* options);

extern "C" void IterativeSolverLinearEquationsInitialize(size_t n, size_t nroot, size_t* range_begin, size_t* range_end,
                                                         const double* rhs, double aughes, double thresh,
                                                         double thresh_value, int hermitian, int verbosity,
                                                         const char* fname, int64_t fcomm, const char* algorithm,
                                                         const char* options);

extern "C" void IterativeSolverNonLinearEquationsInitialize(size_t n, size_t* range_begin, size_t* range_end,
                                                            double thresh, int verbosity, const char* fname,
                                                            int64_t fcomm, const char* algorithm, const char* options);

extern "C" void IterativeSolverOptimizeInitialize(size_t n, size_t* range_begin, size_t* range_end, double thresh,
                                                  double thresh_value, int verbosity, int minimize, const char* fname,
                                                  int64_t fcomm, const char* algorithm, const char* options);

extern "C" void IterativeSolverFinalize();

extern "C" size_t IterativeSolverAddVector(size_t buffer_size, double* parameters, double* action, int sync);

extern "C" void IterativeSolverSolution(int nroot, int* roots, double* parameters, double* action, int sync);

extern "C" size_t IterativeSolverAddValue(double value, double* parameters, double* action, int sync);

extern "C" size_t IterativeSolverEndIteration(size_t buffer_size, double* solution, double* residual, int sync);

extern "C" size_t IterativeSolverAddP(size_t buffer_size, size_t nP, const size_t* offsets, const size_t* indices,
                                      const double* coefficients, const double* pp, double* parameters, double* action,
                                      int sync, void (*func)(const double*, double*, const size_t, const size_t*));

extern "C" void IterativeSolverErrors(double* errors);

extern "C" void IterativeSolverEigenvalues(double* eigenvalues);

extern "C" void IterativeSolverWorkingSetEigenvalues(double* eigenvalues);

extern "C" size_t IterativeSolverSuggestP(const double* solution, const double* residual, size_t maximumNumber,
                                          double threshold, size_t* indices);

extern "C" void IterativeSolverPrintStatistics();

extern "C" int IterativeSolverNonLinear();

extern "C" int IterativeSolverHasValues();
extern "C" int IterativeSolverHasEigenvalues();

extern "C" void IterativeSolverSetDiagonals(const double* diagonals);

extern "C" void IterativeSolverDiagonals(double* diagonals);

extern "C" double IterativeSolverValue();

extern "C" int IterativeSolverVerbosity();

extern "C" int IterativeSolverMaxIter();
extern "C" void IterativeSolverSetMaxIter(int max_iter);

extern "C" int64_t mpicomm_self();

extern "C" int64_t mpicomm_global();
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_ITERATIVESOLVERC_H_
