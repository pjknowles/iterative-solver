!> @example LinearEigensystemExampleF-Pspace-adaptive.F90
!> This is an example of use of the LinearEigensystem framework for iterative
!> finding of the lowest few eigensolutions of a large matrix.
!> A P-space is discovered.
PROGRAM Linear_Eigensystem_Example
  USE Iterative_Solver

  INTEGER, PARAMETER :: n = 1000, nroot = 3, nP = 100
  DOUBLE PRECISION, DIMENSION (n, n) :: m
  DOUBLE PRECISION, DIMENSION (n, nroot) :: c, g
  DOUBLE PRECISION, DIMENSION(nP, nroot) :: p
  DOUBLE PRECISION, DIMENSION (nroot) :: e, error
  LOGICAL, DIMENSION(nroot) :: active
  INTEGER, DIMENSION(0 : nP) :: offsets
  INTEGER, DIMENSION(nP) :: indices
  DOUBLE PRECISION, DIMENSION(nP) :: coefficients
  DOUBLE PRECISION, DIMENSION(nP * nP) :: pp
  INTEGER :: i, j, root
  PRINT *, 'Fortran binding of IterativeSolver'
  m = 1
  DO i = 1, n
    m(i, i) = 3 * i
  END DO

  WRITE (6, *) 'P-space=', nP, ', dimension=', n, ', roots=', nroot
  CALL Iterative_Solver_Linear_Eigensystem_Initialize(n, nroot, thresh = 1d-8, verbosity = 1)
  active = .TRUE.
  CALL Iterative_Solver_Option('convergence', 'residual') ! convergence threshold applies to norm of residual
  offsets(0) = 0
  DO i = 1, nroot
    offsets(i) = i
    indices(i) = i ! the first nroot components
    coefficients(i) = 1
  END DO
  DO i = 1, nroot
    DO j = 1, nroot
      pp(i + (j - 1) * nroot) = m(indices(i), indices(j))
    END DO
  END DO
  CALL Iterative_Solver_Add_P(nroot, offsets, indices, coefficients, pp, c, g, p(: nroot, :))
  g = 0
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
  print *, 'suggest_P returns ', indices(nroot + 1 : nroot + newp)
  DO i = 1, newp
    offsets(nroot + i) = offsets(nroot + i - 1) + 1
    coefficients(nroot + i) = 1
    do j = 1, nroot + newp
      pp(j + (nroot + newp) * (i - 1)) = m(indices(nroot + i), indices(j))
    end do
  END DO
  CALL Iterative_Solver_Add_P(newp, offsets(nroot :) - offsets(nroot), indices(nroot + 1 : nroot + newp), &
      coefficients(nroot + 1 : nroot + newp), pp, c, g, p)
  g = 0
  e = Iterative_Solver_Eigenvalues()
  DO root = 1, nroot
  END DO
  DO iter = 1, 10
    e = Iterative_Solver_Eigenvalues()
    DO root = 1, nroot
      DO i = 1, nP
        DO j = 1, n
          g(j, root) = g(j, root) + m(j, indices(i)) * p(i, root)
        END DO
      END DO
    END DO
    DO root = 1, nroot
      DO j = 1, n
        c(j, root) = c(j, root) - g(j, root) / (m(j, j) - e(i) + 1e-15)
      END DO
    END DO
    IF (Iterative_Solver_End_Iteration(c, g, error, active)) EXIT
    g = MATMUL(m, c)
    CALL Iterative_Solver_Add_Vector(c, g, active, p)
  END DO
  CALL Iterative_Solver_Finalize
  DO i = 1, nroot
    PRINT *, 'solution ', c(1 : MIN(n, 10), i)
  END DO
END PROGRAM Linear_Eigensystem_Example
