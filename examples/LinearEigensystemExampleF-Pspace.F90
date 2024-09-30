module mod_matrix_problem
  use Iterative_Solver, only : mpi_init, mpi_finalize, current_problem
  use iso_c_binding, only : c_int
  use Iterative_Solver_Matrix_Problem, only : matrix_problem
  use Iterative_Solver_Pspace
  INTEGER, PARAMETER :: n = 200, nroot = 3, nP = 30
  DOUBLE PRECISION, DIMENSION (n, n), target :: m
  type(matrix_problem), target :: problem
contains

end module mod_matrix_problem
!> @examples LinearEigensystemExampleF-Pspace-mpi.F90
!> This is an examples of use of the LinearEigensystem framework for iterative
!> finding of the lowest few eigensolutions of a large matrix.
!> A P-space is explicitly declared.
PROGRAM Linear_Eigensystem_Example
  USE Iterative_Solver
  USE iso_c_binding, only : c_funloc
  USE mod_matrix_problem
  DOUBLE PRECISION, DIMENSION (n, nroot) :: c, g
  DOUBLE PRECISION, DIMENSION(nP, nroot) :: p
  DOUBLE PRECISION, DIMENSION (nroot) :: e, error
  DOUBLE PRECISION, DIMENSION(nP) :: coefficients
  DOUBLE PRECISION, DIMENSION(nP, nP) :: pp
  INTEGER :: i, j, root
  LOGICAL :: update
  PRINT *, 'Fortran binding of IterativeSolver'
  CALL MPI_INIT
  call problem%attach(m)
  m = 1
  DO i = 1, n
    m(i, i) = 3 * i
  END DO
  WRITE (6, *) 'P-space=', nP, ', dimension=', n, ', roots=', nroot
  CALL Iterative_Solver_Linear_Eigensystem_Initialize(n, nroot, thresh = 1d-8, verbosity = 0, hermitian = .true.)
!  offsets(0) = 0
!  DO i = 1, nP
!    offsets(i) = i
!    indices(i) = i ! the first nP components
!    coefficients(i) = 1
!  END DO
  !  write (6,*) 'offsets ',offsets
  !  write (6,*) 'indices ',indices
  !  write (6,*) 'coefficients ',coefficients
  call problem%p_space%add_simple([(i, i = 1, nP)]) ! the first nP components
  !  p_space%offsets = offsets
  !  p_space%indices = indices
  !  write (6,*) 'offsets ',p_space%offsets
  !  write (6,*) 'indices ',p_space%indices
  !  write (6,*) 'coefficients ',p_space%coefficients
  !  DO i = 1, problem%p_space%size
  !    DO j = 1, problem%p_space%size
  !      pp(i, j) = m(problem%p_space%indices(i), problem%p_space%indices(j))
  !    END DO
  !  END DO
  !  write (6,*) 'pp',pp
  nwork = nroot
  if (problem%p_space%size .gt.0) then
    current_problem => problem
    nwork = Iterative_Solver_Add_P(nP, problem%p_space%offsets, problem%p_space%indices, problem%p_space%coefficients, problem%pp_action_matrix(), c, g, apply_p_current_problem, .true.)
  else
    c = 0
    do i = 1, nroot
      c(i, i) = 1
    end do
    g = MATMUL(m, c)
    nwork = Iterative_Solver_Add_Vector(c, g)
  end if
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
  CALL MPI_Finalize
CONTAINS
END PROGRAM Linear_Eigensystem_Example
