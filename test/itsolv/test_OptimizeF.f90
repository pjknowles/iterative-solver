! TODO write this
subroutine test_OptimizeF(matrix, n, np, nroot, hermitian, expected_eigenvalues) BIND(C)
  use iso_c_binding
  use Iterative_Solver
  implicit none
  integer(c_size_t), intent(in), value :: n, np, nroot
  integer(c_int), intent(in), value :: hermitian
  double precision, intent(in), dimension(n, n) :: matrix
  double precision, intent(in), dimension(nroot) :: expected_eigenvalues

  double precision, dimension(n, nroot) :: c, g
  double precision, dimension(nroot) :: eigs
  doubleprecision :: eigk, error
  integer :: nwork, i, j, k
  integer, dimension(nroot) :: guess
  double precision :: guess_value

  if (np .gt. 0) return
  write (6, *) 'test_linearEigensystemF ', hermitian
  call Iterative_Solver_Linear_Eigensystem_Initialize(n, nroot, verbosity = 2, hermitian = hermitian.ne.0, &
      thresh = 1d-8)
  nwork = nroot
  c = 0
  do i = 1, nroot
    guess_value = 1d50
    guess(i) = 1
    do k = 1, n
      if (matrix(k, k) .lt. guess_value) then
        do j = 1, i - 1
          if (guess(j).eq.k) goto 44
        end do
        guess(i) = k
        guess_value = matrix(k, k)
      end if
      44 continue
    end do
    c(guess(i), i) = 1
!    write (6, *) 'guess ', guess(i), matrix(guess(i), guess(i))
  end do
  do i = 1, 1000
    g = matmul(matrix, c)
    nwork = Iterative_Solver_Add_Vector(c, g);
    eigs = Iterative_Solver_Working_Set_Eigenvalues(nroot)
!    write (6, *) 'working set eigenvalues ', eigs(:nwork)
    do k = 1, nwork
      eigk = eigs(k)
      do j = 1, n
        g(j, k) = -g(j, k) / (matrix(j, j) + 1e-12 - eigk)
      end do
    end do
    nwork = Iterative_Solver_End_Iteration(c, g);
!    write (6, *) 'error ', Iterative_Solver_Errors()
    if (nwork.le.0) exit
  end do
  eigs = Iterative_Solver_Eigenvalues()
  error = sqrt(dot_product(eigs - expected_eigenvalues, eigs - expected_eigenvalues))
  if (error.gt.1d-10) then
    write (6, *) 'test_linearEigensystemF eigenvalues ', eigs
    write (6, *) 'test_linearEigensystemF expected eigenvalues ', expected_eigenvalues
    write (6, *) 'difference ', eigs - expected_eigenvalues
    write (6, *) 'difference ', sqrt(dot_product(eigs - expected_eigenvalues, eigs - expected_eigenvalues))
    stop
  end if
  call Iterative_Solver_Finalize

end subroutine test_OptimizeF