!> @examples LinearEigensystemExampleF-Pspace-mpi.F90
!> This is an examples of use of the LinearEigensystem framework for iterative
!> finding of the lowest few eigensolutions of a large matrix.
!> A P-space is explicitly declared.

MODULE Pspace
  USE, INTRINSIC :: iso_c_binding
  INTEGER, PARAMETER :: n = 60, nroot = 3, nP = 20
  DOUBLE PRECISION, DIMENSION (n, n) :: m
  INTEGER, DIMENSION(nP) :: indices
  INTEGER :: i, j, root, offset
  CONTAINS
    subroutine apply_on_p(p, g, w_set_size, ranges) BIND(C)
      DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: g
      DOUBLE PRECISION, DIMENSION(nP,nroot), INTENT(inout) :: p
      INTEGER, DIMENSION(*), INTENT(in) :: ranges
      INTEGER(c_size_t), INTENT(in), VALUE :: w_set_size
      INTEGER :: irange, root, range
      !write(*,*) "APPLY_ON_P was called!!!"
      irange = 1
      DO root = 1, nroot ! w_set_size ???
        offset = ranges(irange)
        !range = ranges(irange+1) - ranges(irange)
        !if (rank == 1) then
        !  write(*,*) "gg[",root,"] : ", g((root-1)*n+1:(root-1)*n+range)
        !  write(*,*) "p[",root,"] : ", p(:,root)
        !end if
        DO i = 1, nP
          DO j = ranges(irange)+1, ranges(irange+1)
            !! To be used if action vector is fully stored on each process
            g((root-1)*n+j-offset) = g((root-1)*n+j-offset) + m(j, indices(i)) * p(i,root)
            !! To be used if action vector is distributed across processes
            !g((root-1)*range+j-offset) = g((root-1)*range+j-offset) + m(j, indices(i)) * p(i,root)
          END DO
        END DO
        !if (rank == 1) then
        !  write(*,*) "gg[",root,"] after : ", g((root-1)*n+1:(root-1)*n+range)
        !end if
        irange = irange + 2
      END DO
    end subroutine apply_on_p
END MODULE Pspace

PROGRAM Linear_Eigensystem_Example
  USE Pspace
  USE Iterative_Solver
  USE ProfilerF
  include 'mpif.h'

  !INTEGER, PARAMETER :: n = 20, nroot = 3, nP = 10
  !DOUBLE PRECISION, DIMENSION (n, n) :: m
  DOUBLE PRECISION, DIMENSION (n, nroot) :: c, g
  DOUBLE PRECISION, DIMENSION(nP, nroot) :: p
  DOUBLE PRECISION, DIMENSION (nroot) :: e, error
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: we
  INTEGER, DIMENSION(0 : nP) :: offsets
  !INTEGER, DIMENSION(nP) :: indices
  DOUBLE PRECISION, DIMENSION(nP) :: coefficients
  DOUBLE PRECISION, DIMENSION(nP, nP) :: pp
  !INTEGER :: i, j, root
  LOGICAL :: update
  INTEGER :: nwork, alloc_stat
  INTEGER :: rank, comm_size, ierr
  TYPE(Profiler) :: prof
  rank = 0
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, comm_size, ierr)
  if (rank == 0) then
    PRINT *, 'Fortran binding of IterativeSolver'
    PRINT *, 'Using parallel version'
  endif
  m = 1
  DO i = 1, n
    m(i, i) = 3 * i
  END DO
  prof=Profiler('Eigensystem_Example_P',MPI_COMM_WORLD)
  if (rank == 0) then
    WRITE (6, *) 'P-space=', nP, ', dimension=', n, ', roots=', nroot
  end if
  CALL Iterative_Solver_Linear_Eigensystem_Initialize(n, nroot, pname = 'Eigensystem_Example_P', pcomm = MPI_COMM_WORLD, &
                                                                                  thresh = 1d-8, verbosity = 1)
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
  nwork =  Iterative_Solver_Add_P(nP, offsets, indices, coefficients, pp, c, g, fproc=apply_on_p)
  !g = 0.0d0
  DO iter = 1, 100
    !IF (rank == 0) THEN
    !  PRINT *, 'NWORK:', nwork
    !END IF
    allocate(we(nwork), stat=alloc_stat)
    we = Iterative_Solver_Working_Set_Eigenvalues(nwork)
    !IF (rank == 0) THEN
    !  PRINT *, 'we:', we
    !END IF
    DO root = 1, nwork
      DO j = 1, n
        g(j, root) = - g(j, root) * 1.0d0 / (m(j, j) - we(root) + 1e-15)
      END DO
    END DO
    deallocate(we)
    IF (Iterative_Solver_End_Iteration(c, g, error) == 0) THEN
      EXIT
    END IF
    g = MATMUL(m, c)
    nwork = Iterative_Solver_Pspace_Add_Vector(c, g, fproc=apply_on_p)
  END DO
  CALL Iterative_Solver_Finalize
  call MPI_FINALIZE(ierr)
END PROGRAM Linear_Eigensystem_Example
