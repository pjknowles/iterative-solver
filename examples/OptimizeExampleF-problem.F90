!> @examples OptimizeExampleF-problem.F90
!> This is an example of use of the Optimize framework for iterative
!> minimisation of a non-linear function using the simplified driver.
!> The first example makes stationary a normalised quadratic form, so is equivalent to finding an eigenvector.
!> The second example makes stationary a quadratic form plus a linear force.
module QuasiNewton_Examples
  USE Iterative_Solver_Problem
  private
  !> @brief objective function is (1/2) * c . m . c - sum(c)  where m(i,j) = 1 + (3*i-1)*delta(i,j)
  type, extends(Problem), public :: forced_t
    integer :: size
  contains
    procedure, pass :: residual => forced_residual
    procedure, pass :: diagonals => forced_diagonals
  end type forced_t

  !> @brief objective function is (1/2) * c . m . c / c . c  where m(i,j) = 1 + (3*i-1)*delta(i,j)
  type, extends(Problem), public :: quadratic_t
  contains
    procedure, pass :: residual => quadratic_residual
    procedure, pass :: diagonals => quadratic_diagonals
  end type quadratic_t


contains

  function forced_residual(this, parameters, residuals, range) result(e)
    class(forced_t), intent(in) :: this
    double precision :: e
    double precision, intent(in), dimension(:, :) :: parameters
    double precision, intent(inout), dimension(:, :) :: residuals
    integer, dimension(2), intent(in) :: range
    do i = range(1)+1, range(2); residuals(i, 1) = sum(parameters(:, 1)) + (3 * i - 1) * parameters(i, 1)-1;
    enddo
    e = 0.5 * dot_product(parameters(:, 1), residuals(:, 1)) - 0.5 * sum(parameters(:, 1))
    residuals = residuals - 1
  end function forced_residual

  logical function forced_diagonals(this, d)
    class(forced_t), intent(in) :: this
    double precision, intent(inout), dimension(:) :: d
    d = [(3 * i, i = 1, size(d))]
    forced_diagonals = .true.
  end function forced_diagonals

  function quadratic_residual(this, parameters, residuals, range) result(e)
    class(quadratic_t), intent(in) :: this
    double precision :: e
    double precision, intent(in), dimension(:, :) :: parameters
    double precision, intent(inout), dimension(:, :) :: residuals
    integer, dimension(2), intent(in) :: range
    do i = range(1)+1, range(2); residuals(i, 1) = sum(parameters(:, 1)) + (3 * i - 1) * parameters(i, 1);
    enddo
    e = dot_product(parameters(:, 1), residuals(:, 1)) / dot_product(parameters(:, 1), parameters(:, 1))
    residuals = (residuals - e * parameters) / dot_product(parameters(:, 1), parameters(:, 1))
  end function quadratic_residual

  logical function quadratic_diagonals(this, d)
    class(quadratic_t), intent(in) :: this
    double precision, intent(inout), dimension(:) :: d
    d = [(3 * i, i = 1, size(d))]
    quadratic_diagonals = .true.
  end function quadratic_diagonals

end module QuasiNewton_Examples

PROGRAM QuasiNewton_Example
  USE Iterative_Solver
  USE QuasiNewton_Examples
  IMPLICIT NONE
  INTEGER, PARAMETER :: n = 5, verbosity = 2
  DOUBLE PRECISION, DIMENSION (n) :: c, g
  ! try one of the following
  !  type(quadratic_t) :: problem
  type(forced_t) :: problem

  call mpi_init
  PRINT *, 'Fortran binding of IterativeSolver::IOptimize'

  CALL Iterative_Solver_Optimize_Initialize(n, thresh = 1d-6, verbosity = verbosity, algorithm = "BFGS", &
    options = "max_size_qspace=3")
  c = 0;  c(1) = 1
  CALL Iterative_Solver_Solve(c, g, problem)
  if (verbosity.lt.1) then
    print *, 'Optimized function value ', Iterative_Solver_Value()
    print *, 'Error ', Iterative_Solver_Errors()
  end if
  if (verbosity.gt.1) then
    call Iterative_Solver_Solution([1], c, g)
    PRINT *, 'solution ', c(1:MIN(n, 10))
  end if
  CALL Iterative_Solver_Finalize

  call mpi_finalize
END PROGRAM QuasiNewton_Example
