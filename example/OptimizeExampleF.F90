!> @example OptimizeExampleF.F90
!> This is an example of simplest use of the Optimize framework for iterative
!> solution of non-linear equations
!> The example makes stationary a quadratic form, so is equivalent to finding an eigenvector
PROGRAM QuasiNewton_Example
  USE Iterative_Solver
  IMPLICIT NONE
  INTEGER, PARAMETER :: n = 10
  DOUBLE PRECISION, DIMENSION (n, n) :: m
  DOUBLE PRECISION, DIMENSION (n) :: c, g
  DOUBLE PRECISION :: e, e0
  DOUBLE PRECISION, DIMENSION(1) :: error
  INTEGER :: i, j
  LOGICAL :: converged
  LOGICAL, PARAMETER :: forced = .TRUE.
  PRINT *, 'Fortran binding of IterativeSolver::Optimize'
  m = 1; DO i = 1, n; m(i, i) = i;
  END DO
  CALL Iterative_Solver_Optimize_Initialize(n, thresh = 1d-11, verbosity = 1, algorithm = "BFGS")
  c = 0; if (.not. forced) c(1) = 1
  e0 = m(1, 1)
  DO i = 1, 30
    !    c = c / sqrt(dot_product(c, c))
!    write (6, *) 'c ', c
    g = MATMUL(m, c)
    if (forced) then
      e = 0.5 * dot_product(c, g) - sum(c)
      g = g - 1
    else
      e = dot_product(c, g) / dot_product(c, c)
      g = (g - e * c) / dot_product(c, c)
    end if
    write (6, *) 'e ', e
!    write (6, *) 'g ', g
    CALL Iterative_Solver_Add_Vector(c, g)
!    write (6, *) 'after add_vector c ', c
!    write (6, *) 'after add_vector g ', g
    if (forced) then
      c = c - g / [(m(j, j), j = 1, n)]
    else
      c = c - g / ([(m(j, j), j = 1, n)] - e + 1d-15)
    end if
!    write (6, *) 'c ', c
    converged = Iterative_Solver_End_Iteration(c, g, error)
!    write (6, *) 'after end_iteration c ', c
    IF (converged) EXIT
  END DO
  PRINT *, 'error =', error
  if (.not. forced) c = c / sqrt(dot_product(c, c))
  PRINT *, 'solution ', c(1:MIN(n, 10))
  print *, 'expectation ', e
  CALL Iterative_Solver_Finalize
END PROGRAM QuasiNewton_Example
