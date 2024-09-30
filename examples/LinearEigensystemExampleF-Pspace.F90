!> @examples LinearEigensystemExampleF.F90
!> This is an examples of use of the LinearEigensystem framework for iterative
!> finding of the lowest few eigensolutions of a large matrix.
!> Comparison is made of explicitly declaring a (bad) P-space, generating one automatically, or not using one
PROGRAM Linear_Eigensystem_Example
  USE Iterative_Solver, only : mpi_init, mpi_finalize, mpi_rank_global, Solve_Linear_Eigensystem, Iterative_Solver_Print_Statistics, Iterative_Solver_Finalize
  USE Iterative_Solver_Matrix_Problem, only : matrix_problem
  USE iso_fortran_env, only : output_unit
  INTEGER, PARAMETER :: n = 300, nroot = 3
  INTEGER :: nP, max_p
  DOUBLE PRECISION, DIMENSION (:, :), allocatable, target :: m
  INTEGER :: i
  CALL MPI_INIT
  IF (mpi_rank_global() .gt. 0) close(output_unit)
  ALLOCATE(m(n, n))
  m = 1
  DO i = 1, n
    m(i, i) = 3 * (n + 1 - i)
  END DO
  do nP = 0, 50, 50
    do max_p = 0, 50, 10
      if (max_p.gt.0 .and. nP.gt.0) cycle
      WRITE (6, *) 'Explicit P-space=', nP, ', auto P-space=', max_p, ', dimension=', n, ', roots=', nroot
      call solve(m, nP, max_p)
    end do
  end do
  CALL MPI_Finalize
CONTAINS
  subroutine solve(m, nP, max_P)
    DOUBLE PRECISION, DIMENSION (:, :), INTENT(IN), target :: m
    integer, intent(in) :: nP, max_P
    DOUBLE PRECISION, DIMENSION (n, nroot) :: c, g
    TYPE(Matrix_Problem) :: problem
    logical :: success
    CALL problem%attach(m)
    CALL problem%p_space%add_simple([(i, i = 1, nP)]) ! the first nP components, so not the best
    success = Solve_Linear_Eigensystem(c, g, problem, nroot, verbosity = 2, thresh = 1d-8, hermitian = .true., max_p = max_p)
    if (mpi_rank_global().eq.0) CALL Iterative_Solver_Print_Statistics
    CALL Iterative_Solver_Finalize
  end subroutine solve
END PROGRAM Linear_Eigensystem_Example
