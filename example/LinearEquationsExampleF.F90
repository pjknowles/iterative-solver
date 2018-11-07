!> @example LinearEquationsExampleF.F90
!> This is an example of simplest use of the LinearEquations framework for iterative
!> finding of the lowest few eigensolutions of a large matrix.
PROGRAM Linear_Equations_Example
  USE Iterative_Solver
  IMPLICIT NONE
  INTEGER, PARAMETER :: n = 300, nroot = 2
  DOUBLE PRECISION, PARAMETER :: alpha = 300
  DOUBLE PRECISION, DIMENSION(5), PARAMETER :: augmented_hessian_factors = [0.0_8, .001_8, .01_8, .1_8, 1.0_8]
  DOUBLE PRECISION, DIMENSION (n, n) :: m
  DOUBLE PRECISION, DIMENSION (n, nroot) :: c, g, rhs
  LOGICAL, DIMENSION(nroot) :: active
  DOUBLE PRECISION, DIMENSION (nroot) :: error
  DOUBLE PRECISION :: augmented_hessian
  INTEGER :: i, j, root, iaug
  LOGICAL :: converged
  PRINT *, 'Fortran binding of IterativeSolver'
  DO i = 1, n; m(i, i) = alpha * i + 2 * i - 2; DO j = 1, n; IF (i.NE.j) m(i, j) = i + j - 2;
  END DO;
  END DO
  DO i = 1, nroot; DO j = 1, n; rhs(j, i) = 1 / DBLE(j + i - 1);
  END DO;
  END DO
  DO iaug = 1, SIZE(augmented_hessian_factors)
    augmented_hessian = augmented_hessian_factors(iaug)
    PRINT *, 'solve linear system with augmented hessian factor ', augmented_hessian
    CALL Iterative_Solver_Linear_Equations_Initialize(n, nroot, rhs, augmented_hessian, thresh = 1d-11, verbosity = 1)
    active = .TRUE.
    c = 0; DO i = 1, nroot; c(i, i) = 1;
    ENDDO
    DO i = 1, n
      g = MATMUL(m, c)
      CALL Iterative_Solver_Add_Vector(c, g, active)
      DO root = 1, nroot
        DO j = 1, n
          c(j, root) = c(j, root) - g(j, root) / m(j, j)
        END DO
      END DO
      converged = Iterative_Solver_End_Iteration(c, g, error, active)
      IF (converged) EXIT
    END DO
    PRINT *, 'error =', error
    DO i = 1, nroot
      PRINT *, 'solution ', c(1 : MIN(n, 10), i)
    END DO
    CALL Iterative_Solver_Finalize
  ENDDO
END PROGRAM Linear_Equations_Example
