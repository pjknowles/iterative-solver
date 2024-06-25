!> @examples LinearEigensystemExampleF-problem.F90

program Eigenproblem_Example
  use Iterative_Solver
  use Iterative_Solver_Problem, only: matrix_problem
  implicit none
  double precision, dimension (5000, 5000), target :: matrix=1d0
  double precision, dimension (ubound(matrix,1), 5) :: c, g
  integer :: i, j

  do i = lbound(matrix, 1), ubound(matrix, 1)
    do j = lbound(matrix, 1), ubound(matrix, 1)
      matrix(j, i) = 1
    end do
    matrix(i, i) = 3 * i
  end do

  call Iterative_Solver_Linear_Eigensystem_Initialize(ubound(c, 1), thresh = 1d-7, verbosity = 2, &
    options = "max_size_qspace=10", nroot = ubound(c, 2))
  call Iterative_Solver_Solve(c, g, matrix_problem(matrix), max_iter = 30, generate_initial_guess = .true.)
  print *, 'Eigenvalues ', Iterative_Solver_Eigenvalues()
  if (Iterative_Solver_Verbosity().lt.1) print *, 'Error ', Iterative_Solver_Errors()
  if (Iterative_Solver_Verbosity().gt.1) then
    call Iterative_Solver_Solution([(i, i = lbound(c, 2), ubound(c, 2))], c, g)
    do i = lbound(c, 2), ubound(c, 2)
      print *, 'Eigenvector ', c(1:MIN(ubound(c, 1), 5), i)
    end do
  end if
  call Iterative_Solver_Print_Statistics
  call Iterative_Solver_Finalize

end program Eigenproblem_Example
