!> @brief IterativeSolver Fortran binding
MODULE IterativeSolverF
CONTAINS
 SUBROUTINE IterativeSolverFTest() BIND(C,name='IterativeSolverFTest')
 PRINT *, 'Test Fortran binding of IterativeSolver'
 END SUBROUTINE IterativeSolverFTest
END MODULE IterativeSolverF
