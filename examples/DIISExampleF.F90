!> @example DIISExampleF.F90
!> This is an example of simplest use of the DIIS framework for iterative
!> solution of non-linear equations
!> The example makes stationary a quadratic form, so is equivalent to finding an eigenvector
PROGRAM DIIS_Example
  USE Iterative_Solver
  IMPLICIT NONE
  INTEGER, PARAMETER :: n = 1000
  DOUBLE PRECISION, DIMENSION (n, n) :: m
  DOUBLE PRECISION, DIMENSION (n) :: c, g
  DOUBLE PRECISION :: e
  DOUBLE PRECISION, DIMENSION(1) :: error
  INTEGER :: i, j
  LOGICAL :: converged
  PRINT *, 'Fortran binding of IterativeSolver::DIIS'
  m = 1; DO i = 1, n; m(i, i) = 3 * i;
  END DO
  CALL Iterative_Solver_DIIS_Initialize(n, thresh = 1d-11, verbosity = 1)
  c = 0; c(1) = 1
  DO i = 1, n
    c = c / sqrt(dot_product(c, c))
    g = MATMUL(m, c)
    e = dot_product(c, g)
    g = g - e * c
    CALL Iterative_Solver_Add_Vector(c, g)
    c = c - g / ([(m(j, j), j = 1, n)] - e + 1d-15)
    converged = Iterative_Solver_End_Iteration(c, g, error)
    IF (converged) EXIT
  END DO
  PRINT *, 'error =', error
  c = c / sqrt(dot_product(c, c))
  PRINT *, 'solution ', c(1 : MIN(n, 10))
  print *, 'expectation ', e
  CALL Iterative_Solver_Finalize
END PROGRAM DIIS_Example
