!> @examples LinearEquationsExampleF.F90
!> This is an example of simplest use of the LinearEquations framework for iterative
!> solution of linear equations
module mod_linear_problem
  use Iterative_Solver_Problem, only : Problem
  DOUBLE PRECISION, PARAMETER :: alpha = 300
  type, extends(Problem) :: linear_problem
    integer :: n, nroot
  contains
    procedure, pass :: action
    procedure, pass :: RHS
    procedure, pass :: diagonals
  end type linear_problem
contains
  subroutine action(this, parameters, actions, range)
    class(linear_problem), intent(in) :: this
    double precision, intent(in), dimension(:, :) :: parameters
    double precision, intent(inout), dimension(:, :) :: actions
    integer, dimension(2), intent(in) :: range
    integer :: i, j, k
    do k = lbound(parameters, 2), ubound(parameters, 2)
      do i = range(1) + 1, range(2)
        actions(i, k) = alpha * i * parameters(i, k)
        do j = 1, this%n
          actions(i, k) = actions(i, k) + dble(i + j - 2) * parameters(j, k)
        end do
      end do
    end do
  end subroutine action
  logical function diagonals(this, d)
    class(linear_problem), intent(in) :: this
    double precision, intent(inout), dimension(:) :: d
    integer :: i
    do i = 1, this%n
      d(i) = alpha * i + 2 * i - 2
    end do
    diagonals = .true.
  end function diagonals
  logical function RHS(this, vector, instance, range)
    class(linear_problem), intent(in) :: this
    double precision, intent(inout), dimension(:) :: vector
    integer, intent(in) :: instance
    integer, dimension(2), intent(in) :: range
    integer :: i
    RHS = .false.
    if (instance.lt.1 .or. instance.gt.this%nroot) return
    RHS = .true.
    do i = range(1) + 1, range(2)
      vector(i) = 1 / dble(i + instance - 1)
    end do
  end function RHS
end module mod_linear_problem

PROGRAM Linear_Equations_Example
  USE Iterative_Solver
  USE mod_linear_problem, only : linear_problem
  IMPLICIT NONE
  INTEGER, PARAMETER :: n = 30, nroot = 2
  DOUBLE PRECISION, DIMENSION (n, nroot) :: c, g
  TYPE(linear_problem) :: problem
  call mpi_init
  IF (mpi_rank_global() .gt. 0) close(6)
  PRINT *, 'Fortran binding of IterativeSolver'
  problem%n = n
  problem%nroot = nroot
  CALL Solve_Linear_Equations(c, g, problem, thresh = 1d-11, verbosity = 2)
  PRINT*, 'convergence?', Iterative_Solver_Converged(), ', residual length: ', norm2(g)
  CALL Iterative_Solver_Print_Statistics
  CALL Iterative_Solver_Finalize
  CALL mpi_finalize
END PROGRAM Linear_Equations_Example
