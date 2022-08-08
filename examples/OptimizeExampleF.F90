!> @examples OptimizeExampleF.F90
!> This is an examples of simplest use of the Optimize framework for iterative
!> solution of non-linear equations.
!> With forced=.false., the examples makes stationary a normalised quadratic form, so is equivalent to finding an eigenvector.
!> With forced=.true., the examples makes stationary a quadratic form plus a linear force.
PROGRAM QuasiNewton_Example
  USE Iterative_Solver
  IMPLICIT NONE
  INTEGER, PARAMETER :: n = 2
  DOUBLE PRECISION, DIMENSION (n, n) :: m
  DOUBLE PRECISION, DIMENSION (n) :: c, g
  DOUBLE PRECISION :: e, e0
  INTEGER :: i, j
  LOGICAL :: converged
  LOGICAL, PARAMETER :: forced = .FALSE.
  call mpi_init
  PRINT *, 'Fortran binding of IterativeSolver::IOptimize'
  m = 1; DO i = 1, n; m(i, i) = 3 * i;
  END DO
  CALL Iterative_Solver_Optimize_Initialize(n, thresh = 1d-9, verbosity = 1, algorithm = "BFGS", options="max_size_qspace=7")
  c = 0; if ( .not. forced) c(1) = 1
  e0 = m(1, 1)
  DO i = 1, 30
    g = MATMUL(m, c)
    if (forced) then
      e = 0.5 * dot_product(c, g) - sum(c)
      g = g - 1
    else
      e = dot_product(c, g) / dot_product(c, c)
      g = (g - e * c) / dot_product(c, c)
    end if
    write (6, *) 'function value ', e
    write (6, *) 'c ', c
    write (6, *) 'g ', g
    IF (Iterative_Solver_Add_Vector(c, g, value = e).gt.0) THEN
      if (forced) then
        g = g / [(m(j, j), j = 1, n)]
      else
        g = g / ([(m(j, j), j = 1, n)] - e + 1d-15) &
          + (sum([(c(j) * g(j), j = 1, n)]) / sum([(c(j)**2, j = 1, n)])) * c &
            / ([(m(j, j), j = 1, n)] - e + 1d-15)
      end if
    END IF
    converged = Iterative_Solver_End_Iteration(c, g) .eq.0
    IF (converged) EXIT
  END DO
  if (.not. forced) c = c / sqrt(dot_product(c, c))
  PRINT *, 'solution ', c(1:MIN(n, 10))
  CALL Iterative_Solver_Finalize
  call mpi_finalize
END PROGRAM QuasiNewton_Example
