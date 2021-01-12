!> @examples LinearEigensystemExampleF.F90
!> This is an examples of simplest use of the LinearEigensystem framework for iterative
!> finding of the lowest few eigensolutions of a large matrix.
PROGRAM Linear_Eigensystem_Example
  USE Iterative_Solver
  USE ProfilerF
  include 'mpif.h'
  !USE mpi_f08
  INTEGER, PARAMETER :: n = 60, nroot = 3
  DOUBLE PRECISION, DIMENSION (n, n) :: m
  DOUBLE PRECISION, DIMENSION (n, nroot) :: c, g
  DOUBLE PRECISION, DIMENSION (nroot) :: e
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
  p=Profiler('Eigensystem_Example', 1, 0)
  CALL Iterative_Solver_Linear_Eigensystem_Initialize(n, nroot, thresh = 1d-7, thresh_value = 1d-14, verbosity = 1, &
                                                      pname = 'Eigensystem_Example', mpicomm = MPI_COMM_WORLD)
  c = 0
  DO i = 1, nroot
    c(i, i) = 1;
  END DO
  DO i = 1, n
!    IF (rank == 0) THEN
!      PRINT *, 'ITERATION #', i
!    END IF
    g = MATMUL(m, c)
    nwork = Iterative_Solver_Add_Vector(c, g)
!    IF (rank == 0) THEN
!      PRINT *, 'nwork after Add_Vector():', nwork
!    END IF
    if (nwork == 0) exit
    allocate(we(nwork), stat=alloc_stat)
    we = Iterative_Solver_Working_Set_Eigenvalues(nwork)
!    IF (rank == 0) THEN
!      PRINT *, 'Working set roots after Add_Vector:', we
!    END IF
    DO root = 1, nwork
      DO j = 1, n
        g(j, root) = g(j, root) * 1.0d0 / (m(j, j) - we(root) + 1e-15)
      END DO
    END DO
    deallocate(we)
    nwork = Iterative_Solver_End_Iteration(c, g)
!    IF (rank == 0) THEN
!      PRINT *, 'nwork after End_Iteration():', nwork
!    END IF
    if (nwork == 0) exit
    allocate(we(nwork), stat=alloc_stat)
    we = Iterative_Solver_Working_Set_Eigenvalues(nwork)
!    IF (rank == 0) THEN
!      PRINT *, 'Working set roots after End_Iteration:', we
!    END IF
    deallocate(we)
  END DO
  allocate(we(nroot), stat=alloc_stat)
  we = Iterative_Solver_Eigenvalues()
  IF (rank == 0) THEN
    PRINT *, 'Converged roots:', we
  END IF
  deallocate(we)
  CALL Iterative_Solver_Finalize
  call p%print(6)
  call p%print(6, cumulative = .false.)
  call p%destroy()
  call MPI_FINALIZE(ierr)
END PROGRAM Linear_Eigensystem_Example
