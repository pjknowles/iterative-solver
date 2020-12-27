! Created by Peter Knowles on 26/12/2020.

subroutine test_LinearEquationsF(matrix, rhs, n, np, nroot, hermitian) BIND(C)
  use iso_c_binding
  use Iterative_Solver
  implicit none
  integer(c_size_t), intent(in), value :: n, np, nroot
  integer(c_int), intent(in), value :: hermitian
  double precision, intent(in), dimension(n, n) :: matrix
  double precision, intent(in), dimension(n, nroot) :: rhs

  double precision, dimension(n, nroot) :: c, g
  double precision :: error
  integer :: nwork, i, j, k
  integer, dimension(nroot) :: guess
  double precision :: guess_value

  return ! TODO remove when fortran implementation is complete
  if (np .gt. 0) return
  write (6, *) 'test_linearEquationsF ', hermitian
  call Iterative_Solver_Linear_Equations_Initialize(n, nroot, rhs, hermitian = hermitian.ne.0, &
      thresh = 1d-8)
  nwork = nroot
  c = rhs
  do i = 1, 1000
    g = matmul(matrix, c) - rhs
    nwork = Iterative_Solver_Add_Vector(c, g);
    do k = 1, nwork
      do j = 1, n
        g(j, k) = -g(j, k) / matrix(j, j)
      end do
    end do
    nwork = Iterative_Solver_End_Iteration(c, g);
    !    write (6, *) 'error ', Iterative_Solver_Errors()
    if (nwork.le.0) exit
  end do
  call Iterative_Solver_Solution([(i, i = 1, nroot)], c, g)
  error = 0
  do i = 1, nroot
    error = max(error, sqrt(dot_product(g(:, i), g(:, i))))
    g(:, i) = matmul(matrix, c(:, i)) - rhs(:, i)
    error = max(error, sqrt(dot_product(g(:, i), g(:, i))))
  end do
  if (error.gt.1d-10) then
    write (6, *) 'test_linearEquationsF has failed ', error
    stop
  end if
  call Iterative_Solver_Finalize

end subroutine test_LinearEquationsF