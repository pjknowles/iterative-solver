!> @example OptimizeExampleF.F90
!> This is an example of simplest use of the Optimize framework for iterative
!> solution of non-linear equations
!> The example makes stationary a quadratic form, so is equivalent to finding an eigenvector
PROGRAM QuasiNewton_Example
  USE Iterative_Solver
  IMPLICIT NONE
  INTEGER, PARAMETER :: n = 1000
  DOUBLE PRECISION, DIMENSION (n, n) :: m
  DOUBLE PRECISION, DIMENSION (n) :: c, g
  DOUBLE PRECISION :: e,e0
  DOUBLE PRECISION, DIMENSION(1) :: error
  INTEGER :: i, j
  LOGICAL :: converged
  PRINT *, 'Fortran binding of IterativeSolver::Optimize'
  m = 1; DO i = 1, n; m(i, i) = 3 * i;
  END DO
  CALL Iterative_Solver_Optimize_Initialize(n, thresh = 1d-11, verbosity = 8, algorithm="null")
  c = 0; c(1) = 1
  e0 = m(1,1)
  DO i = 1, n
    c = c / sqrt(dot_product(c, c))
    g = MATMUL(m, c)
    e = dot_product(c, g)
    g = g - e * c
    write (6,*) 'c ',c
    write (6,*) 'g ',g
    CALL Iterative_Solver_Add_Vector(c, g)
    c = c - g / ([(m(j, j), j = 1, n)] - e0 + 1d-15)
    write (6,*) 'c ',c
    converged = Iterative_Solver_End_Iteration(c, g, error)
    IF (converged) EXIT
  END DO
  PRINT *, 'error =', error
  c = c / sqrt(dot_product(c, c))
  PRINT *, 'solution ', c(1 : MIN(n, 10))
  print *, 'expectation ', e
  CALL Iterative_Solver_Finalize
END PROGRAM QuasiNewton_Example
