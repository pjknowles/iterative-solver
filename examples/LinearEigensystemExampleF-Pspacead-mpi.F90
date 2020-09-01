!> @examples LinearEigensystemExampleF-Pspace-adaptive.F90
!> This is an examples of use of the LinearEigensystem framework for iterative
!> finding of the lowest few eigensolutions of a large matrix.
!> A P-space is discovered.
PROGRAM Linear_Eigensystem_Example_P
  USE Iterative_Solver
  USE ProfilerF
  include 'mpif.h'
  INTEGER, PARAMETER :: n = 20, nroot = 3, nP = 10
  DOUBLE PRECISION, DIMENSION (n, n) :: m
  DOUBLE PRECISION, DIMENSION (n, nroot) :: c, g
  DOUBLE PRECISION, DIMENSION (nP, nroot) :: p
  DOUBLE PRECISION, DIMENSION (nroot) :: e, error
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: we
  INTEGER, DIMENSION(0 : nP) :: offsets
  INTEGER, DIMENSION(nP) :: indices
  DOUBLE PRECISION, DIMENSION(nP) :: coefficients
  DOUBLE PRECISION, DIMENSION(nP * nP) :: pp
  INTEGER :: i, j, root, newp
  INTEGER :: nwork, alloc_stat
  INTEGER :: rank, comm_size, ierr
  TYPE(Profiler) :: prof

  call MPI_INIT(ierr)
  !print *, ierr
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  print *, rank
  call MPI_COMM_SIZE(MPI_COMM_WORLD, comm_size, ierr)
  if (rank == 0) then
    PRINT *, 'Fortran binding of IterativeSolver'
    PRINT *, 'Using parallel version'
  endif
  m = 1
  DO i = 1, n
    m(i, i) = 3 * i
  END DO
  prof=Profiler('Eigensystem_Example_P_adapt',MPI_COMM_WORLD)

  !if (rank == 0) then
    WRITE (6, *) 'P-space=', nP, ', dimension=', n, ', roots=', nroot
  !end if
  CALL Iterative_Solver_Linear_Eigensystem_Initialize(n, nroot, pname = 'Eigensystem_Example_P_adapt', &
                                                      pcomm = MPI_COMM_WORLD,thresh = 1d-8, verbosity = 1)
  !CALL Iterative_Solver_Option('convergence', 'residual') ! convergence threshold applies to norm of residual
  indices = 0
  offsets = 0
  coefficients = 0
  !offsets(0) = 0
  DO i = 1, nroot
    offsets(i) = i
    indices(i) = i ! the first nroot components
    coefficients(i) = 1
  END DO
  DO i = 1, nroot
    DO j = 1, nroot
      pp(i + (j - 1) * nroot) = m(indices(i), indices(j))
      !pp(i, j) = m(indices(i), indices(j))
    END DO
  END DO
  nwork = Iterative_Solver_Add_P(nroot, offsets, indices, coefficients, pp, c, g, p(: nroot, :))
  if (nwork .NE. nroot) then
    print *, 'NWORK != NROOT after the first AddP call'
  end if
  g = 0.0d0
  !allocate(we(nwork), stat=alloc_stat)
  e = Iterative_Solver_Eigenvalues()
  DO root = 1, nroot
    DO i = 1, nroot
      DO j = 1, n
        g(j, root) = g(j, root) + m(j, indices(i)) * p(i, root)
      END DO
    END DO
    DO i = 1, nroot
      g(indices(i), root) = 0
    END DO
    DO j = 1, n
      c(j, root) = - g(j, root) / (m(j, j) - e(root) + 1d-10)
    END DO
  END DO
  newp = Iterative_Solver_Suggest_P(c, g, indices(nroot + 1 : nP), 1d-8)
  !if (rank == 0) then
    print *, 'suggest_P returns ', indices(nroot + 1 : nroot + newp), newp
  !end if
  DO i = 1, newp
    offsets(nroot + i) = offsets(nroot + i - 1) + 1
    coefficients(nroot + i) = 1
    do j = 1, nroot + newp
      pp(j + (nroot + newp) * (i - 1)) = m(indices(nroot + i), indices(j))
    end do
  END DO
  nwork = Iterative_Solver_Add_P(newp, offsets(nroot :) - offsets(nroot), indices(nroot + 1 : nroot + newp), &
      coefficients(nroot + 1 : nroot + newp), pp, c, g, p)
  !if (rank == 0) then
    PRINT *, 'NWORK:', nwork
  !end if
  call flush(6)
  !g = 0
  DO iter = 1, 10
    allocate(we(nwork), stat=alloc_stat)
    we = Iterative_Solver_Working_Set_Eigenvalues(nwork)
    DO root = 1, nwork
      DO i = 1, nP
        DO j = 1, n
          g(j, root) = g(j, root) + m(j, indices(i)) * p(i, root)
        END DO
      END DO
    END DO
    DO root = 1, nwork
      DO j = 1, n
        c(j, root) = c(j, root) - g(j, root) / (m(j, j) - we(root) + 1e-15)
      END DO
    END DO
    deallocate(we)
    !IF (Iterative_Solver_End_Iteration(c, g, error)) EXIT
    g = MATMUL(m, c)
    nwork = Iterative_Solver_Add_Vector(c, g, p)
    IF (rank == 0) THEN
      PRINT *, 'NWORK:', nwork
    END IF
    IF (nwork == 0) THEN
      EXIT
    END IF
  END DO
  CALL Iterative_Solver_Finalize
  !DO i = 1, nroot
  !  PRINT *, 'solution ', c(1 : MIN(n, 10), i)
  !END DO
  call MPI_FINALIZE(ierr)
END PROGRAM Linear_Eigensystem_Example_P
