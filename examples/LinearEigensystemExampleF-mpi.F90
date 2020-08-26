!> @examples LinearEigensystemExampleF.F90
!> This is an examples of simplest use of the LinearEigensystem framework for iterative
!> finding of the lowest few eigensolutions of a large matrix.
PROGRAM Linear_Eigensystem_Example
  USE Iterative_Solver
  USE ProfilerF
  include 'mpif.h'
  !USE mpi_f08
  INTEGER, PARAMETER :: n = 6, nroot = 3
  DOUBLE PRECISION, DIMENSION (n, n) :: m
  DOUBLE PRECISION, DIMENSION (n, nroot) :: c, g
  DOUBLE PRECISION, DIMENSION (nroot) :: e, error
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: we
  INTEGER :: i, j, root
  INTEGER :: nwork, alloc_stat
  INTEGER, DIMENSION (nroot) :: roots
  LOGICAL :: converged
  INTEGER :: rank, comm_size, ierr
  TYPE(Profiler) :: p
  !
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, comm_size, ierr)
  if (rank == 0) then
    PRINT *, 'Fortran binding of IterativeSolver'
    PRINT *, 'Using parallel version'
  endif
  !
  m = 1
  DO i = 1, n
    m(i, i) = 3 * i
  END DO
  p=Profiler('Benchmark',MPI_COMM_WORLD)
  CALL Iterative_Solver_Linear_Eigensystem_Initialize(n, nroot, pname = 'Benchmark', pcomm = MPI_COMM_WORLD, thresh = 1d-7, &
                                                       verbosity = 1)
  c = 0
  DO i = 1, nroot
    c(i, i) = 1;
  END DO
  DO i = 1, n
    g = MATMUL(m, c)
    nwork = Iterative_Solver_Add_Vector(c, g)
    IF (nwork .GT. 0) THEN
      allocate(we(nwork), stat=alloc_stat)
      we = Iterative_Solver_Working_Set_Eigenvalues(nwork)
      DO root = 1, nwork
        DO j = 1, n
          c(j, root) = c(j, root) - g(j, root) / (m(j, j) - we(root) + 1e-15)
        END DO
      END DO
    ELSE
      EXIT
    END IF
    !converged = Iterative_Solver_End_Iteration(c, g, error)
    !IF (converged) EXIT
    IF (rank == 0) THEN
      PRINT *, 'NWORK:', nwork
    END IF
    deallocate(we)
  END DO
  e = Iterative_Solver_Eigenvalues()
  IF (rank == 0) THEN
    PRINT *, 'Eigenvalues:', e
  END IF
  !PRINT *, 'error =', error, ' eigenvalue =', e
  CALL Iterative_Solver_Finalize
END PROGRAM Linear_Eigensystem_Example
