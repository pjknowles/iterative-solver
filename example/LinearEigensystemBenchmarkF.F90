PROGRAM Linear_Eigensystem_Benchmark
  USE Iterative_Solver
  USE ProfilerF
  TYPE(Profiler) :: p
  INTEGER, PARAMETER :: n = 2000000, nroot = 3
  !DOUBLE PRECISION, DIMENSION (n, n) :: m
  DOUBLE PRECISION, DIMENSION (n, nroot) :: c, g
  DOUBLE PRECISION, DIMENSION (nroot) :: e, error
  DOUBLE PRECISION :: su
  INTEGER :: i, j, root
  LOGICAL :: converged, update
  PRINT *, 'Fortran binding of IterativeSolver'
  p=Profiler('Benchmark')
  !m = 1; DO i = 1, n; m(i, i) = 3 * i; END DO
  CALL Iterative_Solver_Linear_Eigensystem_Initialize(n, nroot, thresh = 1d-7, verbosity = 1)
  c = 0; DO i = 1, nroot; c(i, i) = 1;
  ENDDO
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
    converged = Iterative_Solver_End_Iteration(c, g, error)
    IF (converged) EXIT
  END DO
  PRINT *, 'error =', error, ' eigenvalue =', e
  CALL Iterative_Solver_Finalize
  call p%print(6)
END PROGRAM Linear_Eigensystem_Benchmark
