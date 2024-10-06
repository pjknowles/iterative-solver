!> @examples LinearEquationsExampleF-matrix.F90
!> This is an examples of simplest use of the LinearEquations framework for iterative
!> solution of linear equations
PROGRAM Linear_Equations_Example
  USE Iterative_Solver, only : mpi_init, mpi_finalize, mpi_rank_global, &
      Solve_Linear_Equations, Iterative_Solver_Print_Statistics, Iterative_Solver_Finalize, Iterative_Solver_Converged
  USE Iterative_Solver_Matrix_Problem, only : Matrix_Problem
  IMPLICIT NONE
  INTEGER, PARAMETER :: n = 300, nroot = 2
  DOUBLE PRECISION, PARAMETER :: alpha = 300
  DOUBLE PRECISION, DIMENSION (n, n), target :: m
  DOUBLE PRECISION, DIMENSION (n, nroot), target :: rhs
  DOUBLE PRECISION, DIMENSION (n, nroot) :: c, g
  TYPE(Matrix_Problem) :: problem
  CALL mpi_init
  PRINT *, 'Fortran binding of IterativeSolver'
  IF (mpi_rank_global() .gt. 0) CLOSE(6)
  CALL initialise_matrices
  CALL problem%attach(m, rhs)
  CALL Solve_Linear_Equations(c, g, problem, thresh = 1d-11, verbosity = 2, max_p = 30, hermitian = .true.)
  PRINT *, 'convergence?', Iterative_Solver_Converged(), ', residual length: ', norm2(g)
  !    print *,c
  !    print *,g
  CALL Iterative_Solver_Print_Statistics
  CALL Iterative_Solver_Finalize
  CALL mpi_finalize
CONTAINS
  SUBROUTINE initialise_matrices
    INTEGER :: i, j
    DO i = 1, n; m(i, i) = alpha * i + 2 * i - 2; DO j = 1, n; IF (i.NE.j) m(i, j) = i + j - 2;
    END DO;
    END DO
    DO i = 1, nroot; DO j = 1, n; rhs(j, i) = 1 / DBLE(j + i - 1);
    END DO;
    END DO
  END SUBROUTINE initialise_matrices
END PROGRAM Linear_Equations_Example
