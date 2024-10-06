module Iterative_Solver_Matrix_Problem
  use Iterative_Solver_Problem, only : Problem
  private
  !> @brief A specialisation of the Problem class for linear problems in which the kernel matrix is stored in full in an existing array
  type, public, extends(Problem) :: Matrix_Problem
    double precision, pointer, dimension(:, :) :: matrix
    double precision, dimension(:, :), pointer :: m_RHS
  contains
    procedure, pass :: attach
    procedure, pass :: RHS
    procedure, pass :: diagonals
    procedure, pass :: action
    procedure, pass :: pp_action_matrix
    procedure, pass :: p_action
  end type Matrix_Problem

contains

  !> @brief Define the problem in terms of an existing kernel matrix, and (for linear equations only), an existing matrix giving the inhomogeneous parts of the equations
  subroutine attach(this, matrix, RHS)
    class(Matrix_Problem), intent(inout) :: this
    double precision, dimension(:, :), target, optional :: matrix
    double precision, dimension(:, :), target, optional :: RHS
    if (present(matrix)) this%matrix => matrix
    if (present(RHS)) this%m_RHS => RHS
  end subroutine attach

  !> @brief Provide the inhomgeneous part of one of the sets of linear equations
  !> @param vector Will contain the requested RHS on exit
  !> @param instance Which equation set for which the RHS should be provided
  !> @param range The range of the space for which actions should be computed. It's OK to provide also the values outside this range (which will happen in a multiprocessing context), but they will be ignored by the solver.
  !> @return Whether the requested instance exists
  logical function RHS(this, vector, instance, range)
    class(Matrix_Problem), intent(in) :: this
    double precision, intent(inout), dimension(:) :: vector
    integer, dimension(2), intent(in) :: range
    integer, intent(in) :: instance
    RHS = .false.
    if (instance.lt.lbound(this%m_RHS, 2).or.instance.gt.ubound(this%m_RHS, 2)) return
    RHS = .true.
    vector(range(1) + 1:range(2)) = this%m_RHS(range(1) + 1:range(2), instance)
  end function RHS

  logical function diagonals(this, d)
    class(Matrix_Problem), intent(in) :: this
    double precision, intent(inout), dimension(:) :: d
    d = [(this%matrix(i, i), i = lbound(this%matrix, 1), ubound(this%matrix, 1))]
    diagonals = .true.
  end function diagonals

  subroutine action(this, parameters, actions, range)
    class(Matrix_Problem), intent(in) :: this
    double precision, intent(in), dimension(:, :) :: parameters
    double precision, intent(inout), dimension(:, :) :: actions
    integer, dimension(2), intent(in) :: range
    actions(range(1) + 1:range(2), :) = matmul(this%matrix(range(1) + 1:range(2), :), parameters)
  end subroutine action

  !> @brief Calculate the representation of the kernel matrix in the P space
  function pp_action_matrix(this) result(matrix)
    class(matrix_problem), intent(in) :: this
    double precision, dimension(:, :), allocatable :: matrix
    allocate(matrix(this%p_space%size, this%p_space%size))
    do i = 1, this%p_space%size
      do j = 1, this%p_space%size
        matrix(i, j) = 0d0
        do ic = this%p_space%offsets(i - 1) + 1, this%p_space%offsets(i)
          do jc = this%p_space%offsets(j - 1) + 1, this%p_space%offsets(j)
            matrix(i, j) = matrix(i, j) + &
                this%matrix(this%p_space%indices(ic), this%p_space%indices(jc)) * &
                    this%p_space%coefficients(ic) * this%p_space%coefficients(jc)
          end do
        end do
      end do
    end do
  end function pp_action_matrix

  !> @brief Calculate the action of the kernel matrix on a set of vectors in the P space
  !> @param p_coefficients The projection of the vectors onto to the P space
  !> @param actions On exit, the computed action has been added to the original contents
  !> @param range The range of the space for which actions should be computed.
  subroutine p_action(this, p_coefficients, actions, range)
    class(matrix_problem), intent(in) :: this
    double precision, dimension(:, :), intent(in) :: p_coefficients
    double precision, dimension(:, :), intent(inout) :: actions
    integer, dimension(2), intent(in) :: range
    do i = lbound(actions, 2), ubound(actions, 2)
      do k = 1, this%p_space%size
        do kc = this%p_space%offsets(k - 1) + 1, this%p_space%offsets(k)
          do j = range(1) + 1, range(2)
            actions(j, i) = actions(j, i) + &
                this%matrix(j, this%p_space%indices(kc)) * this%p_space%coefficients(kc) * &
                    p_coefficients(k, i)
          end do
        end do
      end do
    end do
  end subroutine p_action

end module Iterative_Solver_Matrix_Problem