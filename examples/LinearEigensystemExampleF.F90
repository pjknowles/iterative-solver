!> @example LinearEigensystemExampleF.F90
!> This is an example of simplest use of the LinearEigensystem framework for iterative
!> finding of the lowest few eigensolutions of a large matrix.
PROGRAM Linear_Eigensystem_Example
  USE Iterative_Solver
  INTEGER, PARAMETER :: n = 6, nroot = 3
  DOUBLE PRECISION, DIMENSION (n, n) :: m
  DOUBLE PRECISION, DIMENSION (n, nroot) :: c, g
  DOUBLE PRECISION, DIMENSION (nroot) :: e, error
  LOGICAL, DIMENSION(nroot) :: active
  INTEGER :: i, j, root
  LOGICAL :: converged
  PRINT *, 'Fortran binding of IterativeSolver'
  m = 1; DO i = 1, n; m(i, i) = 3 * i;
  END DO
  CALL Iterative_Solver_Linear_Eigensystem_Initialize(n, nroot, thresh = 1d-7, verbosity = 1)
  active = .true.
  c = 0; DO i = 1, nroot; c(i, i) = 1;
  ENDDO
  DO i = 1, n
    g = MATMUL(m, c)
    CALL Iterative_Solver_Add_Vector(c, g, active, e)
    e = Iterative_Solver_Eigenvalues()
    DO root = 1, nroot
      DO j = 1, n
        c(j, root) = c(j, root) - g(j, root) / (m(j, j) - e(root) + 1e-15)
      END DO
    END DO
    converged = Iterative_Solver_End_Iteration(c, g, error, active)
    IF (converged) EXIT
  END DO
  PRINT *, 'error =', error, ' eigenvalue =', e
  CALL Iterative_Solver_Finalize
END PROGRAM Linear_Eigensystem_Example
