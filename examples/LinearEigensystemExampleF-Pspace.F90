!> @examples LinearEigensystemExampleF-Pspace.F90
!> This is an examples of use of the LinearEigensystem framework for iterative
!> finding of the lowest few eigensolutions of a large matrix.
!> A P-space is explicitly declared.
PROGRAM Linear_Eigensystem_Example
  USE Iterative_Solver

  INTEGER, PARAMETER :: n = 20, nroot = 3, nP = 10
  DOUBLE PRECISION, DIMENSION (n, n) :: m
  DOUBLE PRECISION, DIMENSION (n, nroot) :: c, g
  DOUBLE PRECISION, DIMENSION(nP, nroot) :: p
  DOUBLE PRECISION, DIMENSION (nroot) :: e, error
  INTEGER, DIMENSION(0 : nP) :: offsets
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
  CALL Iterative_Solver_Option('convergence', 'residual') ! convergence threshold applies to norm of residual
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
  CALL Iterative_Solver_Add_P(nP, offsets, indices, coefficients, pp, c, g, p)
  DO iter = 1, 10
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
    IF (Iterative_Solver_End_Iteration(c, g, error)) EXIT
    !write (6,*) 'error=',error
    !write (6,*) 'solution after end_iteration ',c(:,1)
    g = MATMUL(m, c)
    !write (6,*) 'action before add_vector',g(:,1)
    update = Iterative_Solver_Add_Vector(c, g, p)
  END DO
  CALL Iterative_Solver_Finalize
END PROGRAM Linear_Eigensystem_Example
