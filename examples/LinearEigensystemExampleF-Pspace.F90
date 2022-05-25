module matrix_problem
  INTEGER, PARAMETER :: n = 200, nroot = 3, nP = 30
  DOUBLE PRECISION, DIMENSION (n, n) :: m
  INTEGER, DIMENSION(nP) :: indices
  INTEGER, DIMENSION(0:nP) :: offsets
contains
  subroutine apply_p(p, g, nvec, ranges) bind(c)
    use iso_c_binding
    implicit none
    integer(c_size_t), intent(in), value :: nvec
    real(c_double), dimension(n, nvec), intent(inout) :: g
    real(c_double), dimension(nP, nvec), intent(in) :: p
    integer(c_size_t), dimension(2, nvec), intent(in) :: ranges
    integer :: i, j, k
    do i = 1, nvec
      do k = 1, nP
        do j = ranges(1, i) + 1, ranges(2, i)
          g(j, i) = g(j, i) + m(j, indices(k)) * p(indices(k), i)
        end do
      end do
    end do
  end subroutine apply_p
end module matrix_problem
!> @examples LinearEigensystemExampleF-Pspace-mpi.F90
!> This is an examples of use of the LinearEigensystem framework for iterative
!> finding of the lowest few eigensolutions of a large matrix.
!> A P-space is explicitly declared.
PROGRAM Linear_Eigensystem_Example
  USE Iterative_Solver
  USE iso_c_binding, only : c_funloc
  USE matrix_problem
  DOUBLE PRECISION, DIMENSION (n, nroot) :: c, g
  DOUBLE PRECISION, DIMENSION(nP, nroot) :: p
  DOUBLE PRECISION, DIMENSION (nroot) :: e, error
  DOUBLE PRECISION, DIMENSION(nP) :: coefficients
  DOUBLE PRECISION, DIMENSION(nP, nP) :: pp
  INTEGER :: i, j, root
  LOGICAL :: update
  PRINT *, 'Fortran binding of IterativeSolver'
  m = 1
  DO i = 1, n
    m(i, i) = 3 * i
  END DO

  WRITE (6, *) 'P-space=', nP, ', dimension=', n, ', roots=', nroot
  CALL Iterative_Solver_Linear_Eigensystem_Initialize(n, nroot, thresh = 1d-8, verbosity = 0, hermitian = .true.)
  offsets(0) = 0
  DO i = 1, nP
    offsets(i) = i
    indices(i) = i ! the first nP components
    coefficients(i) = 1
  END DO
  DO i = 1, nP
    DO j = 1, nP
      pp(i, j) = m(indices(i), indices(j))
    END DO
  END DO
  nwork = Iterative_Solver_Add_P(nP, offsets, indices, coefficients, pp, c, g, apply_p, .true.)
  DO iter = 1, n
    e = Iterative_Solver_Eigenvalues()
    write (6, *) 'eigenvalues=', Iterative_Solver_Eigenvalues()
    DO root = 1, nroot
      DO j = 1, n
        c(j, root) = c(j, root) - g(j, root) / (m(j, j) - e(i) + 1e-15)
      END DO
    END DO
    nwork = Iterative_Solver_End_Iteration(c, g)
    write (6, *) 'error=', Iterative_Solver_Errors()
    IF (nwork.le.0) EXIT
    g = MATMUL(m, c)
    nwork = Iterative_Solver_Add_Vector(c, g)
  END DO
  call Iterative_Solver_Print_Statistics
  !  CALL Iterative_Solver_Solution([(i,i=1,nroot)],c,g)
  !  write (6,*) 'final solution ',c
  !  write (6,*) 'final residual ',g
  CALL Iterative_Solver_Finalize
CONTAINS
END PROGRAM Linear_Eigensystem_Example
