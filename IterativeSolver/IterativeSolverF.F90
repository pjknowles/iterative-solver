!> @brief IterativeSolver Fortran binding
MODULE IterativeSolverF
 TYPE Davidson
 END TYPE Davidson
CONTAINS
 SUBROUTINE IterativeSolverFTest() BIND(C,name='IterativeSolverFTest')
 PRINT *, 'Test Fortran binding of IterativeSolver'
 END SUBROUTINE IterativeSolverFTest
 subroutine anharmonic_residual()
 end subroutine anharmonic_residual
END MODULE IterativeSolverF
