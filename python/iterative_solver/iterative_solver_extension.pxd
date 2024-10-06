cdef extern from "molpro/linalg/IterativeSolverC.h":
    void IterativeSolverLinearEigensystemInitialize(size_t nQ, size_t nroot, size_t* range_begin,
                                                    size_t* range_end, double thresh, double thresh_value,
                                                    int hermitian, int verbosity, const char* fname,
                                                    int fcomm, const char* algorithm, const char* options)

    void IterativeSolverLinearEquationsInitialize(size_t n, size_t nroot, size_t* range_begin, size_t* range_end,
                                                  const double* rhs, double aughes, double thresh,
                                                  double thresh_value, int hermitian, int verbosity,
                                                  const char* fname, int fcomm, const char* algorithm,
                                                  const char* options)

    void IterativeSolverNonLinearEquationsInitialize(size_t n, size_t* range_begin, size_t* range_end,
                                                     double thresh, int verbosity, const char* fname,
                                                     int fcomm, const
                                                     char* algorithm, const
                                                     char* options
                                                     )
    void IterativeSolverOptimizeInitialize(size_t n, size_t* range_begin, size_t* range_end, double thresh,
                                           double thresh_value, int verbosity, int minimize, const char* fname,
                                           int fcomm, const char* algorithm, const char* options)

    void IterativeSolverFinalize()

    size_t IterativeSolverAddVector(size_t buffer_size, double* parameters, double* action, int sync)

    void IterativeSolverSolution(int nroot, int* roots, double* parameters, double* action, int sync)

    size_t IterativeSolverAddValue(double value, double* parameters, double* action, int sync)

    size_t IterativeSolverEndIteration(size_t buffer_size, double* solution, double* residual, int sync)

    int IterativeSolverEndIterationNeeded();

    # size_t IterativeSolverAddP(size_t buffer_size, size_t nP, const size_t* offsets, const size_t* indices,
    # const double* coefficients, const double* pp, double* parameters, double* action,
    # int sync, void (*func)(const double*, double*, const size_t, const size_t*))

    void IterativeSolverErrors(double* errors)

    void IterativeSolverEigenvalues(double* eigenvalues)

    void IterativeSolverWorkingSetEigenvalues(double* eigenvalues)

    size_t IterativeSolverSuggestP(const double* solution, const double* residual, size_t maximumNumber,
    double threshold, size_t* indices)

    void IterativeSolverPrintStatistics()

    int IterativeSolverNonLinear()

    int IterativeSolverHasValues()
    int IterativeSolverHasEigenvalues()

    void IterativeSolverSetDiagonals(const double* diagonals)

    void IterativeSolverDiagonals(double* diagonals)

    double IterativeSolverValue()

    int IterativeSolverVerbosity()

    int IterativeSolverMaxIter()
    void IterativeSolverSetMaxIter(int max_iter)

    int IterativeSolver_mpicomm_self()
    int IterativeSolver_mpicomm_global()

    int IterativeSolverConverged()