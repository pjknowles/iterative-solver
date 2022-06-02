module Iterative_Solver_Problem

  private

  !> @brief Abstract class defining the problem-specific interface for the simplified solver
  !> interface to IterativeSolver
  !type, abstract, public :: Problem
  type, public :: Problem
  contains
    procedure, pass :: diagonals
    procedure, pass :: precondition
    procedure, pass :: residual
    procedure, pass :: action
    procedure, pass :: report
  end type Problem

contains

  !> @brief Optionally provide the diagonal elements of the underlying kernel. If
  !> implemented and returning true, the provided diagonals will be used by
  !> IterativeSolver for preconditioning (and therefore the precondition() function does
  !> not need to be implemented), and, in the case of linear problems, for selection of
  !> the P space. Otherwise, preconditioning will be done with precondition(), and any P
  !> space has to be provided manually.
  !> @param d On exit, contains the diagonal elements
  !> @return Whether diagonals have been provided.
  logical function diagonals(this, d)
    class(Problem), intent(in) :: this
    double precision, intent(inout), dimension(:) :: d
    diagonals = .false.
  end function diagonals

  !> @brief Apply preconditioning to a residual vector in order to predict a step towards
  !> the solution
  !> @param residual On entry, assumed to be the residual. On exit, the negative of the
  !> predicted step.
  !> @param shift When called from LinearEigensystem, contains the corresponding current
  !> eigenvalue estimates for each of the parameter vectors in the set. All other solvers
  !> should pass a vector of zeroes, which is the default if omitted.
  !> @param diagonals The diagonal elements of the underlying kernel. If passed, they will be used,
  !> otherwise the default preconditioner does nothing.
  subroutine precondition(this, action, shift, diagonals)
    class(Problem), intent(in) :: this
    double precision, intent(inout), dimension(:, :) :: action
    double precision, intent(in), dimension(:), optional :: shift
    double precision, intent(in), dimension(:), optional :: diagonals
    double precision, parameter :: small = 1e-14
    if (present(diagonals)) then
      do i = lbound(action, 2), ubound(action, 2)
        if (present(shift)) then
          do j = lbound(action, 1), ubound(action, 1)
            action(j, i) = action(j, i) / (diagonals(j) + shift(i) + small)
          end do
        else
          do j = lbound(action, 1), ubound(action, 1)
            action(j, i) = action(j, i) / (diagonals(j) + small)
          end do
        end if
      end do
    end if
  end subroutine precondition


  !> @brief Calculate the residual vector. Used by non-linear solvers (NonLinearEquations,
  !> Optimize) only.
  !> @param parameters The trial solution for which the residual is to be calculated
  !> @param resid The residual vector
  !> @return In the case where the residual is an exact differential, the corresponding
  !> function value. Used by Optimize but not NonLinearEquations.
  function residual(this, parameters, residuals) result(value)
    class(Problem), intent(in) :: this
    double precision :: value
    double precision, intent(in), dimension(:, :) :: parameters
    double precision, intent(inout), dimension(:, :) :: residuals
    residuals = 0d0
  end function residual

  !> @brief Calculate the action of the kernel matrix on a set of parameters. Used by
  !> linear solvers, but not by the non-linear solvers (NonLinearEquations, Optimize).
  !> @param parameters The trial solutions for which the action is to be calculated
  !> @param act The action vectors
  subroutine action(this, parameters, actions)
    class(Problem), intent(in) :: this
    double precision, intent(in), dimension(:, :) :: parameters
    double precision, intent(inout), dimension(:, :) :: actions
  end subroutine action

  !> @brief Report progress at the end of each iteration, or at the end of the calculation
  !> @return .true. if the information was used, and therefore the caller should be silent
  logical function report(this, iteration, verbosity, errors, value, eigenvalues)
    class(Problem), intent(in) :: this
    integer, intent(in) :: iteration !< The iteration number if positive, or zero indicating successful convergence, or a negative number indicating failure to converge
    integer, intent(in) :: verbosity !< Expected interpretation:
    !< - 0 or less, nothing should be printed
    !< - 1 nothing should be printed if iteration > 0
    !< - 2 or more show progress in each iteration
    double precision, intent(in), dimension(:) :: errors !< The current residual norm for each solution
    double precision, intent(in), optional :: value !< In the case of optimisation only, the current objective function value
    double precision, dimension(:), intent(in), optional :: eigenvalues !< In the case of eigenproblem only, the current eigenvalues
    if ((iteration.le.0 .and. verbosity.ge.1) .or. verbosity.ge.2) then
      if (iteration.gt. 0 .and. verbosity.ge.2) then
        write (6, '(A,I3,1X,A,(T32,10F7.2))') 'Iteration', iteration, 'log10(|residual|)=', log10(errors)
      else if (iteration.eq.0) then
        write (6, '(A,(T32,10F7.2))') 'Converged,   log10(|residual|)=', log10(errors)
      else
        write (6, '(A,(T32,10F7.2))') 'Unconverged, log10(|residual|)=', log10(errors)
      end if
      if (present(value)) then
        write (6, *) 'Objective function value ', value
      end if
      if (present(eigenvalues)) then
        write (6, *) 'Eigenvalues ', eigenvalues
      end if
    end if
    report = .true.
  end function report

end module Iterative_Solver_Problem

module try_Iterative_Solver_Problem
  use Iterative_Solver_Problem
  type, extends(Problem) :: my_Problem

  end type my_Problem

contains
  subroutine try
    type(my_Problem) :: thing
  end subroutine try
end module try_Iterative_Solver_Problem