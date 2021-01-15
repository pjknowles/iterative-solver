! TODO write this
! f = 0.5 * (x-b).h.(x-b) where b=[111...]
! df/dx = h x - h.b
function test_OptimizeF(matrix, n) BIND(C)
  use iso_c_binding
  use Iterative_Solver
  implicit none
  integer(c_int) :: test_OptimizeF
  integer(c_size_t), intent(in), value :: n
  double precision, intent(in), dimension(n, n) :: matrix

  double precision, dimension(n) :: c, g
  doubleprecision :: eigk, error, value
  double precision, parameter :: thresh=1d-8
  integer :: nwork, i, j, k

  test_OptimizeF = 1
  write (6, *) 'test_OptimizeF '
  call Iterative_Solver_Optimize_Initialize(n, verbosity = 2, &
      thresh = thresh, mpicomm = mpicomm_compute())
  nwork = 1
  c = 0
  c(1) = 1
  do i = 1, 1000
    g = matmul(matrix, c - 1)
    value = 0.5d0 * dot_product(g, c - 1)
    if(Iterative_Solver_Add_Value(value, c, g)) then
      do j = 1, n
        g(j) = g(j) / matrix(j, j)
      end do
    end if
    nwork = Iterative_Solver_End_Iteration(c, g);
    if (nwork.le.0) exit
  end do
  error = sqrt(dot_product(c - 1, c - 1))
!  write (6,*) 'difference error ',error, thresh
  if (error.gt.thresh) then
    write (6, *) 'OptimizeF has failed ',value
    test_OptimizeF = 0
  end if
  call Iterative_Solver_Finalize

end function test_OptimizeF