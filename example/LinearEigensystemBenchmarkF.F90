PROGRAM Linear_Eigensystem_Benchmark
  USE Iterative_Solver
  USE ProfilerF
#ifdef HAVE_MPI_H
  include 'mpif.h'
#endif
  TYPE(Profiler) :: p
  INTEGER, PARAMETER :: n = 2000000, nroot = 3
  !DOUBLE PRECISION, DIMENSION (n, n) :: m
  DOUBLE PRECISION, DIMENSION (n, nroot) :: c, g
  DOUBLE PRECISION, DIMENSION (nroot) :: e, error
  DOUBLE PRECISION :: su
  INTEGER :: i, j, root, ierr
  LOGICAL :: converged, update
  PRINT *, 'Fortran binding of IterativeSolver'
#ifdef HAVE_MPI_H
  call MPI_Init(ierr)
  PRINT *, 'Using parallel version'
  p=Profiler('Benchmark',MPI_COMM_WORLD)
  CALL Iterative_Solver_Linear_Eigensystem_Initialize(n, nroot, pname = 'Benchmark', pcomm = MPI_COMM_WORLD, thresh = 1d-7, &
                                                      verbosity = 1)
#else
  PRINT *, 'Using serial version'
  p=Profiler('Benchmark')   
  CALL Iterative_Solver_Linear_Eigensystem_Initialize(n, nroot, pname = 'Benchmark', thresh = 1d-7, verbosity = 1)
#endif
  !m = 1; DO i = 1, n; m(i, i) = 3 * i; END DO
  c = 0; DO i = 1, nroot; c(i, i) = 1; ENDDO
  DO i = 1, n
!    g = MATMUL(m, c)
    call p%start('residual')
    do root=1,nroot
      su=0d0
      do j=1,n
    su = su + c(j,root)
      end do
    do j=1,n; g(j,root) = (3*j-1)*c(j,root)+su; end do
      end do
    call p%stop('residual')
    call p%start('Add_Vector')
    update = Iterative_Solver_Add_Vector(c, g, e)
    call p%stop('Add_Vector')
    IF (update) THEN
      call p%start('update')
      e = Iterative_Solver_Eigenvalues()
      DO root = 1, nroot
        DO j = 1, n
          c(j, root) = c(j, root) - g(j, root) / (3*j - e(root) + 1e-15)
        END DO
      END DO
      call p%stop('update')
    END IF
    call p%start('End_Iter')
    converged = Iterative_Solver_End_Iteration(c, g, error)
    call p%stop('End_Iter')
    IF (converged) EXIT
  END DO
  PRINT *, 'error =', error, ' eigenvalue =', e
  call p%start('Finalize')
  CALL Iterative_Solver_Finalize
  call p%stop('Finalize')
  call p%print(6)
#ifdef HAVE_MPI_H
  call MPI_Finalize(ierr)
#endif
END PROGRAM Linear_Eigensystem_Benchmark
