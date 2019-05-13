!> @example OptimizeExampleF.F90
!> This is an example of simplest use of the Optimize framework for iterative
!> solution of non-linear equations.
!> With forced=.false., the example makes stationary a normalised quadratic form, so is equivalent to finding an eigenvector.
!> With forced=.true., the example makes stationary a quadratic form plus a linear force.
PROGRAM QuasiNewton_Example
  USE Iterative_Solver
  IMPLICIT NONE
  INTEGER, PARAMETER :: n = 100
  DOUBLE PRECISION, DIMENSION (n, n) :: m
  DOUBLE PRECISION, DIMENSION (n) :: c, g
  DOUBLE PRECISION :: e, e0
  DOUBLE PRECISION, DIMENSION(1) :: error
  INTEGER :: i, j
  LOGICAL :: converged
  LOGICAL, PARAMETER :: forced = .FALSE.
  PRINT *, 'Fortran binding of IterativeSolver::Optimize'
  m = 1; DO i = 1, n; m(i, i) = 3*i;
  END DO
  CALL Iterative_Solver_Optimize_Initialize(n, thresh = 1d-9, verbosity = 1, algorithm = "L-BFGS")
  c = 0; if (.not. forced) c(1) = 1
  e0 = m(1, 1)
  DO i = 1, 30
        c = c / sqrt(dot_product(c, c))
    g = MATMUL(m, c)
    if (forced) then
      e = 0.5 * dot_product(c, g) - sum(c)
      g = g - 1
    else
      e = dot_product(c, g) / dot_product(c, c)
      g = (g - e * c) / dot_product(c, c)
    end if
    write (6, *) 'function value ',e
    IF (Iterative_Solver_Add_Value(e, c, g)) THEN
      if (forced) then
        c = c - g / [(m(j, j), j = 1, n)]
      else
        c = c - g / ([(m(j, j), j = 1, n)] - e + 1d-15) &
            + (sum([(c(j)*g(j),j=1,n)])/sum([(c(j)**2,j=1,n)])) * c &
                / ([(m(j, j), j = 1, n)] - e + 1d-15)
      end if
    END IF
    converged = Iterative_Solver_End_Iteration(c, g, error)
    IF (converged) EXIT
  END DO
  if (.not. forced) c = c / sqrt(dot_product(c, c))
  PRINT *, 'solution ', c(1:MIN(n, 10))
  CALL Iterative_Solver_Finalize
END PROGRAM QuasiNewton_Example
