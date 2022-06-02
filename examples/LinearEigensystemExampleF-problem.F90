!> @examples LinearEigensystemExampleF-problem.F90
module Example_Problem
  use Iterative_Solver_Problem
  implicit none
  private
  !> @brief matrix m(i,j) = 1 + (diagonal_factor * i - 1) * delta(i,j)
  type, extends(Problem), public :: problem_t
    double precision :: diagonal_factor = 1d0
  contains
    procedure, pass :: action, diagonals
  end type problem_t

contains

  subroutine action(this, parameters, actions)
    class(problem_t), intent(in) :: this
    double precision, intent(in), dimension(:, :) :: parameters
    double precision, intent(inout), dimension(:, :) :: actions
    integer :: i, j
    do j = lbound(actions, 2), ubound(actions, 2)
      do i = lbound(actions, 1), ubound(actions, 1)
        actions(i, j) = sum(parameters(:, j)) + (this%diagonal_factor * i - 1) * parameters(i, j)
      enddo
    enddo
  end subroutine action

  logical function diagonals(this, d)
    class(problem_t), intent(in) :: this
    double precision, intent(inout), dimension(:) :: d
    integer :: i
    d = [(this%diagonal_factor * i, i = 1, size(d))]
    diagonals = .true.
  end function diagonals

end module Example_Problem

program Eigenproblem_Example
  use Iterative_Solver
  use Example_Problem
  implicit none
  interface
    subroutine mpi_init() bind (C, name = 'mpi_init')
    end subroutine mpi_init
    subroutine mpi_finalize() bind (C, name = 'mpi_finalize')
    end subroutine mpi_finalize
  end interface
  double precision, dimension (5000, 10) :: c, g
  integer :: i

  call mpi_init

  call Iterative_Solver_Linear_Eigensystem_Initialize(ubound(c, 1), thresh = 1d-7, verbosity = 2, &
    options = "max_size_qspace=10", nroot = ubound(c, 2))
  call Iterative_Solver_Solve(c, g, problem_t(3d0), max_iter = 30, generate_initial_guess = .true.)
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

  call mpi_finalize
end program Eigenproblem_Example
