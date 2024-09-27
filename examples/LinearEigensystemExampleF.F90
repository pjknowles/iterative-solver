!> @examples LinearEigensystemExampleF.F90
!> This is an example of simplest use of the LinearEigensystem framework for iterative
!> finding of the lowest few eigensolutions of a large matrix.
program Eigenproblem_Example
  use Iterative_Solver, only : Solve_Linear_Eigensystem, Iterative_Solver_Verbosity, Iterative_Solver_Solution, Iterative_Solver_Eigenvalues, Iterative_Solver_Errors, Iterative_Solver_Finalize, Iterative_Solver_Print_Statistics, mpi_init, mpi_rank_global, mpi_finalize
  use Iterative_Solver_Problem, only : matrix_problem
  use iso_fortran_env, only : output_unit
  implicit none
  double precision, dimension (500, 500), target :: matrix = 1d0
  double precision, dimension (ubound(matrix, 1), 5) :: c, g
  integer :: i, j
  logical :: converged
  call mpi_init
  IF (mpi_rank_global() .gt. 0) close(output_unit)
  do i = lbound(matrix, 1), ubound(matrix, 1)
    matrix(i, i) = 3 * i
  end do
  converged = Solve_Linear_Eigensystem(c, g, matrix_problem(matrix) &
      , nroot = ubound(c, 2) &
      , thresh = 1d-7 &
      , verbosity = 2 &
      , max_iter = 30&
      )
  print *, 'Converged?', converged, ', Eigenvalues ', Iterative_Solver_Eigenvalues()
  if (Iterative_Solver_Verbosity().lt.1) print *, 'Error ', Iterative_Solver_Errors()
  if (Iterative_Solver_Verbosity().gt.1) then
    call Iterative_Solver_Solution([(i, i = lbound(c, 2), ubound(c, 2))], c, g)
    do i = lbound(c, 2), ubound(c, 2)
      print *, 'Eigenvector ', c(1:MIN(ubound(c, 1), 5), i)
    end do
  end if
  call Iterative_Solver_Print_Statistics
  call Iterative_Solver_Finalize
  call mpi_finalize
end program Eigenproblem_Example
