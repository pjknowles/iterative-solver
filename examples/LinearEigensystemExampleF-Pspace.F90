!> @examples LinearEigensystemExampleF-Pspace-mpi.F90
!> This is an examples of use of the LinearEigensystem framework for iterative
!> finding of the lowest few eigensolutions of a large matrix.
!> A P-space is explicitly declared.
PROGRAM Linear_Eigensystem_Example
  USE Iterative_Solver
  interface
    subroutine mpi_init() BIND (C, name = 'mpi_init')
    end subroutine mpi_init
    subroutine mpi_finalize() BIND (C, name = 'mpi_finalize')
    end subroutine mpi_finalize
    !    function mpi_comm_global() BIND (C, name = 'mpi_comm_global')
    !      use iso_c_binding, only: c_int64_t
    !      integer(c_int64_t) mpi_comm_global
    !    end function mpi_comm_global
  end interface

  INTEGER, PARAMETER :: n = 20, nroot = 3, nP = 10
  DOUBLE PRECISION, DIMENSION (n, n) :: m
  DOUBLE PRECISION, DIMENSION (n, nroot) :: c, g
  DOUBLE PRECISION, DIMENSION(nP, nroot) :: p
  DOUBLE PRECISION, DIMENSION (nroot) :: e, error
  INTEGER, DIMENSION(0:nP) :: offsets
  INTEGER, DIMENSION(nP) :: indices
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
  CALL Iterative_Solver_Linear_Eigensystem_Initialize(n, nroot, thresh = 1d-8, verbosity = 1)
!  CALL Iterative_Solver_Option('convergence', 'residual') ! convergence threshold applies to norm of residual
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
  nwork = Iterative_Solver_Add_P(nP, offsets, indices, coefficients, pp, c, g, apply_p)
  DO iter = 1, n
    e = Iterative_Solver_Eigenvalues()
    DO root = 1, nroot
      DO i = 1, nP
        DO j = 1, n
          g(j, root) = g(j, root) + m(j, indices(i)) * p(i, root)
        END DO
      END DO
    END DO
    !write (6,*) 'residual after adding p-space contribution ',g(:,1)
    DO root = 1, nroot
      DO j = 1, n
        c(j, root) = c(j, root) - g(j, root) / (m(j, j) - e(i) + 1e-15)
      END DO
    END DO
    !write (6,*) 'solution after update ',c(:,1)
    nwork = Iterative_Solver_End_Iteration(c, g)
    IF (nwork.le.0) EXIT
    !write (6,*) 'error=',error
    !write (6,*) 'solution after end_iteration ',c(:,1)
    g = MATMUL(m, c)
    !write (6,*) 'action before add_vector',g(:,1)
    nwork = Iterative_Solver_Add_Vector(c, g)
  END DO
  CALL Iterative_Solver_Finalize
  CALL mpi_finalize
CONTAINS
  subroutine apply_p(p, g, upd_size, ranges) bind(c)
    use iso_c_binding
    implicit none
    integer(c_size_t), intent(in), value :: upd_size
    real(c_double), dimension(*), intent(inout) :: g
    real(c_double), dimension(nP, upd_size), intent(in) :: p
    integer(c_size_t), dimension(2,upd_size), intent(in) :: ranges
    !  double precision, pointer, contiguous :: temp(:,:) => NULL ()
    integer :: i, j, k
    do i = 1, upd_size
      do k = 1, nP
        do j = ranges(1,i)+1,ranges(2,i)
          g(j) = g(j) + m(j, indices(k)) * p(indices(k),i)
        end do
      end do
    end do

    !  temp => memory_allocate(stsym_data%n_ci,upd_size)
    !  temp = 0.0d0
    !  call mu_casci_Hc_p(stsym_data,zint,p,temp,upd_size)
    !  irange = 1
    !  do i = 1, upd_size
    !    range = ranges(irange+1) - ranges(irange)
    !    g_from = (i-1)*stsym_data%n_ci + 1
    !    g_to = (i-1)*stsym_data%n_ci + range
    !    p_from = ranges(irange) + 1
    !    p_to = ranges(irange+1)
    !    call daxpy_X(range,1d0,temp(p_from:p_to,i),1,g(g_from:g_to),1)
    !    irange = irange + 2
    !  enddo
    !  call memory_release(temp)
  end subroutine apply_p
END PROGRAM Linear_Eigensystem_Example
