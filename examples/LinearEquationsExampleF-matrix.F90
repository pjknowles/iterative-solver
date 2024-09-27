!> @examples LinearEquationsExampleF.F90
!> This is an examples of simplest use of the LinearEquations framework for iterative
!> solution of linear equations
module mod_linear_problem
  use Iterative_Solver_Problem
  type, extends(matrix_Problem) :: linear_problem
    double precision, dimension(:, :), pointer :: rhss
  contains
    procedure, pass :: RHS
  end type linear_problem
contains
  logical function RHS(this, vector, instance)
    class(linear_problem), intent(in) :: this
    double precision, intent(inout), dimension(:) :: vector
    integer, intent(in) :: instance
    RHS = .false.
    if (instance.lt.lbound(this%rhss, 2).or.instance.gt.ubound(this%rhss, 2)) return
    RHS = .true.
    vector = this%rhss(:, instance)
    print *,vector
    print *,(this%matrix(i,i), i=lbound(this%matrix,1), ubound(this%matrix,1))
  end function RHS
end module mod_linear_problem

PROGRAM Linear_Equations_Example
  USE Iterative_Solver
  USE mod_linear_problem, only : linear_problem
  IMPLICIT NONE
  INTEGER, PARAMETER :: n = 30, nroot = 2
  DOUBLE PRECISION, PARAMETER :: alpha = 300
  DOUBLE PRECISION, DIMENSION(1), PARAMETER :: augmented_hessian_factors = [0.0_8]! issue 510 , .001_8, .01_8, .1_8, 1.0_8]
  DOUBLE PRECISION, DIMENSION (n, n), target :: m
  DOUBLE PRECISION, DIMENSION (n, nroot), target :: rhs
  DOUBLE PRECISION, DIMENSION (n, nroot) :: c, g
  DOUBLE PRECISION :: augmented_hessian
  INTEGER :: iaug
  LOGICAL :: converged
  type(linear_problem) :: problem
  call mpi_init
  PRINT *, 'Fortran binding of IterativeSolver'
  if (mpi_rank_global() .gt. 0) close(6)
  call initialise_matrices
  problem = linear_problem(m, rhs)
  DO iaug = 1, SIZE(augmented_hessian_factors)
    augmented_hessian = augmented_hessian_factors(iaug)
    PRINT *, 'solve linear system with augmented hessian factor ', augmented_hessian
    converged = Solve_Linear_Equations(c, g, problem, augmented_hessian = augmented_hessian, thresh = 1d-11, verbosity = 2)
    print*, 'convergence?', converged, ', residual length: ', norm2(g)
    !    print *,c
    !    print *,g
    Call Iterative_Solver_Print_Statistics
    CALL Iterative_Solver_Finalize
  ENDDO
  CALL mpi_finalize
contains
  subroutine initialise_matrices
    INTEGER :: i, j
    DO i = 1, n; m(i, i) = alpha * i + 2 * i - 2; DO j = 1, n; IF (i.NE.j) m(i, j) = i + j - 2; END DO; END DO
    DO i = 1, nroot; DO j = 1, n; rhs(j, i) = 1 / DBLE(j + i - 1); END DO; END DO
  end subroutine initialise_matrices
END PROGRAM Linear_Equations_Example
