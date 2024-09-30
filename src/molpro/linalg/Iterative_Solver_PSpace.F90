module Iterative_Solver_Pspace
  private
  type, public :: PSpace
    integer, dimension(:), allocatable :: indices
    integer, dimension(:), allocatable :: offsets
    double precision, dimension(:), allocatable :: coefficients
    integer :: size = 0 !< the dimension of the P space
    logical :: simple !< whether the space consists entirely of elements that are simple elements of the full space, rather than linear combinations
  contains
    procedure, pass :: add_complex => PSpace_add_complex
    procedure, pass :: add_simple => PSpace_add_simple
    procedure, pass :: ensure => PSpace_ensure
    final :: PSpace_final
  end type PSpace

contains

  subroutine PSpace_ensure(this)
    class(PSpace), intent(inout) :: this
    if (.not. allocated(this%indices)) then
      this%size = 0
      allocate(this%offsets(0:0), this%indices(0), this%coefficients(0))
      this%offsets(0) = 0
      this%simple = .true.
    end if
  end subroutine PSpace_ensure
  subroutine PSpace_add_complex(this, indices, coefficients)
    implicit none
    class(PSpace), intent(inout) :: this
    integer, dimension(:), intent(in) :: indices
    double precision, dimension(:), intent(in) :: coefficients
    integer, dimension(:), allocatable :: tempi
    double precision, dimension(:), allocatable :: tempd
    integer :: n_coeff
    call this%ensure()
    n_coeff = size(this%coefficients)
    this%simple = this%simple .and. n_coeff.le.1
    allocate(tempi(0:this%size))
    tempi(0:this%size) = this%offsets(0:this%size)
    deallocate(this%offsets)
    allocate(this%offsets(0:this%size + 1))
    this%offsets(0:this%size) = tempi
    this%size = this%size + 1
    this%offsets(this%size) = n_coeff + size(coefficients)
    deallocate(tempi)
    allocate(tempi(n_coeff))
    allocate(tempd(n_coeff))
    tempi = this%indices(:n_coeff)
    tempd = this%coefficients(:n_coeff)
    deallocate(this%indices, this%coefficients)
    allocate(this%indices(n_coeff + size(coefficients)))
    allocate(this%coefficients(n_coeff + size(coefficients)))
    this%indices(:n_coeff) = tempi
    this%coefficients(:n_coeff) = tempd
    this%indices(n_coeff + 1:n_coeff + size(coefficients)) = indices
    this%coefficients(n_coeff + 1:n_coeff + size(coefficients)) = coefficients
    deallocate(tempi, tempd)
  end subroutine PSpace_add_complex
  subroutine PSpace_add_simple(this, indices)
    class(PSpace), intent(inout) :: this
    integer, dimension(:), intent(in) :: indices
    do i = lbound(indices, 1), ubound(indices, 1)
      call this%add_complex([indices(i)], [1d0])
    end do
  end subroutine PSpace_add_simple
  subroutine PSpace_final(this)
    type(PSpace) :: this
    if (allocated(this%indices)) deallocate(this%indices, this%offsets, this%coefficients)
  end subroutine PSpace_final

end module Iterative_Solver_Pspace
