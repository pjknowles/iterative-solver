!> @examples LinearEigensystemExampleF-Pspace-mpi.F90
!> This is an examples of use of the LinearEigensystem framework for iterative
!> finding of the lowest few eigensolutions of a large matrix.
!> A P-space is explicitly declared.
!MODULE ABSTRACT
!  ABSTRACT INTERFACE
!    subroutine func_tmpl(a, b) BIND(C)
!      !DOUBLE PRECISION, DIMENSION(*), INTENT(in) :: a, b
!      DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: a, b
!    end subroutine func_tmpl
!  END INTERFACE
!END MODULE

MODULE Pspace
  INTEGER, PARAMETER :: n = 20, nroot = 3, nP = 10
  DOUBLE PRECISION, DIMENSION (n, n) :: m
  INTEGER, DIMENSION(nP) :: indices
  INTEGER :: i, j, root
  CONTAINS
    subroutine apply_on_p(p, g) BIND(C)
      DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: p, g
      !DOUBLE PRECISION, INTENT(in) :: p(:,:)
      !DOUBLE PRECISION, INTENT(in) :: g(:,:)
      write(*,*) "APPLY_ON_P was called!!!"
      !DO root = 1, nwork
      !  DO i = 1, nP
      !    DO j = 1, n
      !      g(j, root) = g(j, root) + m(j, indices(i)) * p(i, root)
      !    END DO
      !  END DO
      !END DO
    end subroutine apply_on_p
END MODULE Pspace

PROGRAM Linear_Eigensystem_Example
  USE Pspace
  !USE ABSTRACT
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
  !PROCEDURE(func_tmpl), POINTER :: func_app => apply_on_p
  TYPE(Profiler) :: prof
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
  !nwork =  Iterative_Solver_Add_P(nP, offsets, indices, coefficients, pp, c, g, p, fproc_ptr=func_app)
  nwork =  Iterative_Solver_Add_P(nP, offsets, indices, coefficients, pp, c, g, p, fproc=apply_on_p)
  g = 0.0d0
  DO iter = 1, 100
    allocate(we(nwork), stat=alloc_stat)
    we = Iterative_Solver_Working_Set_Eigenvalues(nwork)
    DO root = 1, nwork
      DO i = 1, nP
        DO j = 1, n
          g(j, root) = g(j, root) + m(j, indices(i)) * p(i, root)
        END DO
      END DO
    END DO
    !write (6,*) 'residual after adding p-space contribution ',g(:,1)
    DO root = 1, nwork
      DO j = 1, n
        g(j, root) = - g(j, root) * 1.0d0 / (m(j, j) - we(root) + 1e-15)
      END DO
    END DO
    deallocate(we)
    !write (6,*) 'solution after update ',c(:,1)
    !IF (Iterative_Solver_End_Iteration(c, g, error)) EXIT
    !write (6,*) 'error=',error
    !write (6,*) 'solution after end_iteration ',c(:,1)
    g = MATMUL(m, c)
    !write (6,*) 'action before add_vector',g(:,1)
    nwork = Iterative_Solver_Add_Vector(c, g, p)
    IF (rank == 0) THEN
      PRINT *, 'NWORK:', nwork
    END IF
    IF (nwork == 0) THEN
      EXIT
    END IF
  END DO
  CALL Iterative_Solver_Finalize
  call MPI_FINALIZE(ierr)
END PROGRAM Linear_Eigensystem_Example
