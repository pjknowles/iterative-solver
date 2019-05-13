!> @brief IterativeSolver Fortran binding
MODULE Iterative_Solver
  USE iso_c_binding
  PUBLIC :: Iterative_Solver_Linear_Eigensystem_Initialize, Iterative_Solver_Finalize
  PUBLIC :: Iterative_Solver_DIIS_Initialize, Iterative_Solver_Linear_Equations_Initialize
  PUBLIC :: Iterative_Solver_Optimize_Initialize
  PUBLIC :: Iterative_Solver_Add_Vector, Iterative_Solver_End_Iteration
  PUBLIC :: Iterative_Solver_Add_P, Iterative_Solver_Suggest_P
  PUBLIC :: Iterative_Solver_Eigenvalues
  PUBLIC :: Iterative_Solver_Option
  PRIVATE
  INTEGER(c_size_t) :: m_nq, m_nroot

CONTAINS

  !> \brief Finds the lowest eigensolutions of a matrix using Davidson's method, i.e. preconditioned Lanczos.
  !> Example of simplest use: @include LinearEigensystemExampleF.F90
  !> Example including use of P space: @include LinearEigensystemExampleF-Pspace.F90
  SUBROUTINE Iterative_Solver_Linear_Eigensystem_Initialize(nq, nroot, thresh, maxIterations, verbosity, orthogonalize)
    INTEGER, INTENT(in) :: nq !< dimension of matrix
    INTEGER, INTENT(in) :: nroot !< number of eigensolutions desired
    DOUBLE PRECISION, INTENT(in), OPTIONAL :: thresh !< convergence threshold
    INTEGER, INTENT(in), OPTIONAL :: maxIterations !< maximum number of iterations
    INTEGER, INTENT(in), OPTIONAL :: verbosity !< how much to print. Default is zero, which prints nothing except errors.
    !< One gives a single progress-report line each iteration.
    LOGICAL, INTENT(in), OPTIONAL :: orthogonalize !< whether to orthogonalize expansion vectors (default true)
    INTERFACE
      SUBROUTINE Iterative_Solver_Linear_Eigensystem_InitializeC(nq, nroot, thresh, maxIterations, verbosity, orthogonalize) &
          BIND(C, name = 'IterativeSolverLinearEigensystemInitialize')
        USE iso_c_binding
        INTEGER(C_size_t), INTENT(in), VALUE :: nq
        INTEGER(C_size_t), INTENT(in), VALUE :: nroot
        REAL(c_double), INTENT(in), VALUE :: thresh
        INTEGER(C_int), INTENT(in), VALUE :: maxIterations
        INTEGER(C_int), INTENT(in), VALUE :: verbosity
        INTEGER(C_int), INTENT(in), VALUE :: orthogonalize
      END SUBROUTINE Iterative_Solver_Linear_Eigensystem_InitializeC
    END INTERFACE
    INTEGER(c_int) :: verbosityC = 0, maxIterationsC = 0, orthogonalizeC = 1
    REAL(c_double) :: threshC = 0d0
    m_nq = INT(nq, kind = c_size_t)
    m_nroot = INT(nroot, kind = c_size_t)
    IF (PRESENT(thresh)) THEN
      threshC = thresh
    END IF
    IF (PRESENT(maxIterations)) THEN
      maxIterationsC = INT(maxIterations, kind = c_int)
    END IF
    IF (PRESENT(verbosity)) THEN
      verbosityC = INT(verbosity, kind = c_int)
    END IF
    IF (PRESENT(orthogonalize)) THEN
      IF (orthogonalize) orthogonalizeC = 1
      IF (.NOT.orthogonalize) orthogonalizeC = 0
    END IF
    CALL Iterative_Solver_Linear_Eigensystem_InitializeC(m_nq, m_nroot, threshC, maxIterationsC, verbosityC, orthogonalizeC)
  END SUBROUTINE Iterative_Solver_Linear_Eigensystem_Initialize

  !> \brief Finds the solutions of linear equation systems using a generalisation of Davidson's method, i.e. preconditioned Lanczos
  !> Example of simplest use: @include LinearEquationsExampleF.F90
  !! Example including use of P space: include LinearEquationsExampleF-Pspace.F90
  SUBROUTINE Iterative_Solver_Linear_Equations_Initialize(nq, nroot, rhs, augmented_hessian, thresh, maxIterations, &
      verbosity, orthogonalize)
    INTEGER, INTENT(in) :: nq !< dimension of matrix
    INTEGER, INTENT(in) :: nroot !< number of eigensolutions desired
    double precision, INTENT(in), DIMENSION(nq, nroot) :: rhs !< the constant right-hand-side of each equation system
    double precision, INTENT(in), OPTIONAL :: augmented_hessian
    !< If zero, solve the inhomogeneous equations unmodified. If 1, solve instead
    !< the augmented hessian problem. Other values scale the augmented hessian damping.
    DOUBLE PRECISION, INTENT(in), OPTIONAL :: thresh !< convergence threshold
    INTEGER, INTENT(in), OPTIONAL :: maxIterations !< maximum number of iterations
    INTEGER, INTENT(in), OPTIONAL :: verbosity !< how much to print. Default is zero, which prints nothing except errors.
    !< One gives a single progress-report line each iteration.
    LOGICAL, INTENT(in), OPTIONAL :: orthogonalize !< whether to orthogonalize expansion vectors (default true)
    INTERFACE
      SUBROUTINE Iterative_Solver_Linear_Equations_InitializeC(nq, nroot, rhs, augmented_hessian, thresh, maxIterations, &
          verbosity, orthogonalize) &
          BIND(C, name = 'IterativeSolverLinearEquationsInitialize')
        USE iso_c_binding
        INTEGER(C_size_t), INTENT(in), VALUE :: nq
        INTEGER(C_size_t), INTENT(in), VALUE :: nroot
        REAL(c_double), INTENT(in), DIMENSION(*) :: rhs
        REAL(c_double), INTENT(in), VALUE :: augmented_hessian
        REAL(c_double), INTENT(in), VALUE :: thresh
        INTEGER(C_int), INTENT(in), VALUE :: maxIterations
        INTEGER(C_int), INTENT(in), VALUE :: verbosity
        INTEGER(C_int), INTENT(in), VALUE :: orthogonalize
      END SUBROUTINE Iterative_Solver_Linear_Equations_InitializeC
    END INTERFACE
    INTEGER(c_int) :: verbosityC = 0, maxIterationsC = 0, orthogonalizeC = 1
    REAL(c_double) :: threshC = 0d0, augmented_hessianC = 0d0
    m_nq = INT(nq, kind = c_size_t)
    m_nroot = INT(nroot, kind = c_size_t)
    IF (PRESENT(augmented_hessian)) THEN
      augmented_hessianC = augmented_hessian
    END IF
    IF (PRESENT(thresh)) THEN
      threshC = thresh
    END IF
    IF (PRESENT(maxIterations)) THEN
      maxIterationsC = INT(maxIterations, kind = c_int)
    END IF
    IF (PRESENT(verbosity)) THEN
      verbosityC = INT(verbosity, kind = c_int)
    END IF
    IF (PRESENT(orthogonalize)) THEN
      IF (orthogonalize) orthogonalizeC = 1
      IF (.NOT.orthogonalize) orthogonalizeC = 0
    END IF
    CALL Iterative_Solver_Linear_Equations_InitializeC(m_nq, m_nroot, rhs, augmented_hessianC, threshC, maxIterationsC, &
        verbosityC, orthogonalizeC)
  END SUBROUTINE Iterative_Solver_Linear_Equations_Initialize

  !> \brief Optimization
  !> through the L-BFGS or related methods.
  !> Example of simplest use: @include OptimizeExampleF.F90
  SUBROUTINE Iterative_Solver_Optimize_Initialize(nq, thresh, maxIterations, verbosity, algorithm, minimize)
    INTEGER, INTENT(in) :: nq !< dimension of parameter space
    DOUBLE PRECISION, INTENT(in), OPTIONAL :: thresh !< convergence threshold
    INTEGER, INTENT(in), OPTIONAL :: maxIterations !< maximum number of iterations
    INTEGER, INTENT(in), OPTIONAL :: verbosity !< how much to print. Default is zero, which prints nothing except errors.
    CHARACTER(*), INTENT(in), OPTIONAL :: algorithm !< keyword specifying optimization algorithm
    LOGICAL, INTENT(in), OPTIONAL :: minimize !< whether to minimize (default) or maximize
    !< One gives a single progress-report line each iteration.
    INTERFACE
      SUBROUTINE Iterative_Solver_Optimize_InitializeC(nq, thresh, maxIterations, verbosity, algorithm, minimize) &
          BIND(C, name = 'IterativeSolverOptimizeInitialize')
        USE iso_c_binding
        INTEGER(C_size_t), INTENT(in), VALUE :: nq
        REAL(c_double), INTENT(in), VALUE :: thresh
        INTEGER(C_int), INTENT(in), VALUE :: maxIterations
        INTEGER(C_int), INTENT(in), VALUE :: verbosity
        INTEGER(C_int), INTENT(in), VALUE :: minimize
        CHARACTER(kind = c_char), DIMENSION(*), INTENT(in) :: algorithm
      END SUBROUTINE Iterative_Solver_Optimize_InitializeC
    END INTERFACE
    INTEGER(c_int) :: verbosityC, maxIterationsC, minimizeC
    REAL(c_double) :: threshC
    CHARACTER(LEN=128) :: algorith
    m_nq = INT(nq, kind = c_size_t)
    m_nroot = 1
    IF (PRESENT(thresh)) THEN
      threshC = thresh
    ELSE
      threshC = 0
    END IF
    IF (PRESENT(maxIterations)) THEN
      maxIterationsC = INT(maxIterations, kind = c_int)
    ELSE
      maxIterationsC = 0
    END IF
    IF (PRESENT(verbosity)) THEN
      verbosityC = INT(verbosity, kind = c_int)
    ELSE
      verbosityC = 0
    END IF
    IF (PRESENT(algorithm)) THEN
      algorith = algorithm
    ELSE
      algorith = ""
    END IF
    minimize_C = 1
    IF (PRESENT(minimize)) THEN
      IF (.NOT. minimize) minimizeC = 0
    END IF
    CALL Iterative_Solver_Optimize_InitializeC(m_nq, threshC, maxIterationsC, verbosityC, c_string_c(algorith),minimizeC)
  CONTAINS
    FUNCTION c_string_c(fstring)
      CHARACTER(*), INTENT(in) :: fstring
      CHARACTER(kind = c_char), DIMENSION(:), ALLOCATABLE :: c_string_c
      INTEGER :: i
      ALLOCATE(CHARACTER(kind = c_char) :: c_string_c(len_TRIM(fstring) + 1))
      DO i = 1, len_TRIM(fstring)
        c_string_c(i) = fstring(i:i)
      END DO
      c_string_c(len_TRIM(fstring) + 1) = c_null_char
    END FUNCTION c_string_c
  END SUBROUTINE Iterative_Solver_Optimize_Initialize


  !> \brief Accelerated convergence of non-linear equations
  !> through the DIIS or related methods.
  !> Example of simplest use: @include DIISExampleF.F90
  SUBROUTINE Iterative_Solver_DIIS_Initialize(nq, thresh, maxIterations, verbosity)
    INTEGER, INTENT(in) :: nq !< dimension of parameter space
    DOUBLE PRECISION, INTENT(in), OPTIONAL :: thresh !< convergence threshold
    INTEGER, INTENT(in), OPTIONAL :: maxIterations !< maximum number of iterations
    INTEGER, INTENT(in), OPTIONAL :: verbosity !< how much to print. Default is zero, which prints nothing except errors.
    !< One gives a single progress-report line each iteration.
    INTERFACE
      SUBROUTINE Iterative_Solver_DIIS_InitializeC(nq, thresh, maxIterations, verbosity) &
          BIND(C, name = 'IterativeSolverDIISInitialize')
        USE iso_c_binding
        INTEGER(C_size_t), INTENT(in), VALUE :: nq
        REAL(c_double), INTENT(in), VALUE :: thresh
        INTEGER(C_int), INTENT(in), VALUE :: maxIterations
        INTEGER(C_int), INTENT(in), VALUE :: verbosity
      END SUBROUTINE Iterative_Solver_DIIS_InitializeC
    END INTERFACE
    INTEGER(c_int) :: verbosityC, maxIterationsC
    REAL(c_double) :: threshC
    m_nq = INT(nq, kind = c_size_t)
    m_nroot = 1
    IF (PRESENT(thresh)) THEN
      threshC = thresh
    ELSE
      threshC = 0
    END IF
    IF (PRESENT(maxIterations)) THEN
      maxIterationsC = INT(maxIterations, kind = c_int)
    ELSE
      maxIterationsC = 0
    END IF
    IF (PRESENT(verbosity)) THEN
      verbosityC = INT(verbosity, kind = c_int)
    ELSE
      verbosityC = 0
    END IF
    CALL Iterative_Solver_DIIS_InitializeC(m_nq, threshC, maxIterationsC, verbosityC)
  END SUBROUTINE Iterative_Solver_DIIS_Initialize

  !> \brief Terminate the iterative solver
  SUBROUTINE Iterative_Solver_Finalize
    INTERFACE
      SUBROUTINE IterativeSolverFinalize() BIND(C, name = 'IterativeSolverFinalize')
        USE iso_c_binding
      END SUBROUTINE IterativeSolverFinalize
    END INTERFACE
    CALL IterativeSolverFinalize
  END SUBROUTINE Iterative_Solver_Finalize

  !> \brief Take, typically, a current solution and residual, add it to the expansion set, and return new solution.
  !> In the context of Lanczos-like linear methods, the input will be a current expansion vector and the result of
  !> acting on it with the matrix, and the output will be a new expansion vector.
  !> \param parameters On input, the current solution or expansion vector. On exit, the interpolated solution vector.
  !> \param action On input, the residual for parameters (non-linear), or action of matrix on parameters (linear).
  !> On exit, the expected (non-linear) or actual (linear) residual of the interpolated parameters.
  !> \param parametersP On exit, the interpolated solution projected onto the P space.
  !> \return whether it is expected that the client should make an update, based on the returned parameters and residual, before the subsequent call to Iterative_Solver_End_Iteration()
  FUNCTION Iterative_Solver_Add_Vector(parameters, action, parametersP)
    USE iso_c_binding
    LOGICAL :: Iterative_Solver_Add_Vector
    DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: parameters
    DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: action
    DOUBLE PRECISION, DIMENSION(*), INTENT(inout), OPTIONAL :: parametersP
    INTERFACE
      FUNCTION Iterative_Solver_Add_Vector_C(parameters, action, parametersP) &
          BIND(C, name = 'IterativeSolverAddVector')
        USE iso_c_binding
        INTEGER(c_int) Iterative_Solver_Add_Vector_C
        REAL(c_double), DIMENSION(*), INTENT(inout) :: parameters
        REAL(c_double), DIMENSION(*), INTENT(inout) :: action
        REAL(c_double), DIMENSION(*), INTENT(inout) :: parametersP
      END FUNCTION Iterative_Solver_Add_Vector_C
    END INTERFACE
    DOUBLE PRECISION, DIMENSION(0) :: pdummy
    IF (PRESENT(parametersP)) THEN
      Iterative_Solver_Add_Vector = Iterative_Solver_Add_Vector_C(parameters, action, parametersP).NE.0
    ELSE
      Iterative_Solver_Add_Vector = Iterative_Solver_Add_Vector_C(parameters, action, pdummy).NE.0
    END IF
  END FUNCTION Iterative_Solver_Add_Vector

  !>@brief Take the updated solution vector set, and adjust it if necessary so that it becomes the vector to
  !> be used in the next iteration; this is done only in the case of linear solvers where the orthogonalize option is set.
  !> Also calculate the degree of convergence, and write progress to standard output
  !> \param solution The current solution, after interpolation and updating with the preconditioned residual.
  !> \param residual The residual after interpolation.
  !> \param error Error indicator for each sought root.
  !> \return .TRUE. if convergence reached for all roots
  LOGICAL FUNCTION Iterative_Solver_End_Iteration(solution, residual, error)
    USE iso_c_binding
    DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: solution
    DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: residual
    DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: error
    INTERFACE
      INTEGER(c_int) FUNCTION Iterative_Solver_End_Iteration_C(solution, residual, error) &
          BIND(C, name = 'IterativeSolverEndIteration')
        USE iso_c_binding
        REAL(c_double), DIMENSION(*), INTENT(inout) :: solution
        REAL(c_double), DIMENSION(*), INTENT(inout) :: residual
        REAL(c_double), DIMENSION(*), INTENT(inout) :: error
      END FUNCTION Iterative_Solver_End_Iteration_C
    END INTERFACE
    Iterative_Solver_End_Iteration = &
        Iterative_Solver_End_Iteration_C(solution, residual, error) /= 0
  END FUNCTION Iterative_Solver_End_Iteration


  !> \brief add P-space vectors to the expansion set, and return new solution.
  !> \param nP the number of P-space vectors to add
  !> \param offsets specifies the start point in indices and coefficients that defines each vector.
  !> \param indices Index in the full space of a contribution to a new P vector
  !> \param coefficients Value of a contribution to a new P vector
  !> \param pp The P-P block of the matrix, dimensioned (number of existing P + nP, nP)
  !> \param parameters On input, the current solution or expansion vector. On exit, the interpolated solution vector.
  !> \param action On input, the residual for parameters (non-linear), or action of matrix on parameters (linear).
  !> On exit, the expected (non-linear) or actual (linear) residual of the interpolated parameters.
  !> \param parametersP On exit, the interpolated solution projected onto the P space.
  SUBROUTINE Iterative_Solver_Add_P(nP, offsets, indices, coefficients, pp, parameters, action, parametersP)
    INTEGER, INTENT(in) :: nP
    INTEGER, INTENT(in), DIMENSION(0 : nP) :: offsets
    INTEGER, INTENT(in), DIMENSION(offsets(nP)) :: indices
    DOUBLE PRECISION, DIMENSION(offsets(nP)), INTENT(in) :: coefficients
    DOUBLE PRECISION, DIMENSION(*), INTENT(in) :: pp
    DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: parameters
    DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: action
    DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: parametersP
    INTERFACE
      SUBROUTINE IterativeSolverAddPC(nP, offsets, indices, coefficients, pp, parameters, action, parametersP) &
          BIND(C, name = 'IterativeSolverAddP')
        USE iso_c_binding
        INTEGER(c_size_t), INTENT(in), VALUE :: nP
        INTEGER(c_size_t), INTENT(in), DIMENSION(0 : nP) :: offsets
        INTEGER(c_size_t), INTENT(in), DIMENSION(offsets(nP)) :: indices
        REAL(c_double), DIMENSION(offsets(nP)), INTENT(in) :: coefficients
        REAL(c_double), DIMENSION(*), INTENT(in) :: pp
        REAL(c_double), DIMENSION(*), INTENT(inout) :: parameters
        REAL(c_double), DIMENSION(*), INTENT(inout) :: action
        REAL(c_double), DIMENSION(*), INTENT(inout) :: parametersP
      END SUBROUTINE IterativeSolverAddPC
    END INTERFACE
    INTEGER(c_size_t), DIMENSION(0 : nP) :: offsetsC
    INTEGER(c_size_t), DIMENSION(SIZE(indices)) :: indicesC
    offsetsC = INT(offsets, c_size_t)
    !write (6,*) 'fortrann addp nP ',nP
    !write (6,*) 'fortrann addp offsets ',offsets
    !write (6,*) 'fortrann addp offsetsC ',offsetsC
    !write (6,*) 'indices ',indices
    do i = 1, offsets(nP)
      indicesC(i) = INT(indices(i) - 1, c_size_t) ! 1-base to 0-base
      !write (6,*) 'fortran addp',indicesC(i)
    end do
    !write (6,*) 'indicesC ',indicesC
    CALL IterativeSolverAddPC(INT(nP, c_size_t), offsetsC, indicesC, coefficients, &
        pp, parameters, action, parametersP)
  END SUBROUTINE Iterative_Solver_Add_P

  !> \brief Take an existing solution and its residual, and suggest P vectors
  !> \param solution On input, the current solution.
  !> \param residual On input, the residual for solution.
  !> \param indices On exit, the most important base vectors
  !> \param threshold Base vectors whose predicted contribution is less than
  !> than this are not considered
  !> \return The number of vectors suggested.
  FUNCTION Iterative_Solver_Suggest_P(solution, residual, indices, threshold)
    INTEGER :: Iterative_Solver_Suggest_P
    DOUBLE PRECISION, DIMENSION(*), INTENT(in) :: solution
    DOUBLE PRECISION, DIMENSION(*), INTENT(in) :: residual
    INTEGER, INTENT(inout), DIMENSION(:) :: indices
    DOUBLE PRECISION, INTENT(in), OPTIONAL :: threshold
    INTERFACE
      FUNCTION IterativeSolverSuggestP(solution, residual, maximumNumber, threshold, indices) &
          BIND(C, name = 'IterativeSolverSuggestP')
        USE iso_c_binding
        INTEGER(c_size_t) :: IterativeSolverSuggestP
        INTEGER(c_size_t), VALUE :: maximumNumber
        INTEGER(c_size_t), INTENT(inout), DIMENSION(maximumNumber) :: indices
        REAL(c_double), DIMENSION(*), INTENT(in) :: solution
        REAL(c_double), DIMENSION(*), INTENT(in) :: residual
        REAL(c_double), INTENT(in), VALUE :: threshold
      END FUNCTION IterativeSolverSuggestP
    END INTERFACE
    REAL(C_double) :: thresholdC = 0
    INTEGER(c_size_t), DIMENSION(SIZE(indices)) :: indicesC
    INTEGER(c_size_t) :: maximumNumber
    maximumNumber = INT(size(indices), c_size_t)
    !write (6,*) 'fortran suggestP, maximumNumber=',size(indices)
    IF (PRESENT(threshold)) thresholdC = threshold
    Iterative_Solver_Suggest_P = INT(&
        IterativeSolverSuggestP(solution, residual, maximumNumber, thresholdC, indicesC) &
        )
    do i = 1, Iterative_Solver_Suggest_P
      indices(i) = int(indicesC(i)) + 1
    end do
  END FUNCTION Iterative_Solver_Suggest_P

  !> \brief give options to the iterative solver
  SUBROUTINE Iterative_Solver_Option(key, val)
    CHARACTER(len = *), INTENT(in) :: key !< Option name
    CHARACTER(len = *), INTENT(in) :: val !< Option value
    INTERFACE
      SUBROUTINE IterativeSolverOption(key, val) BIND(C, name = 'IterativeSolverOption')
        USE iso_c_binding
        CHARACTER(kind = c_char), DIMENSION(*), INTENT(in) :: key, val
      END SUBROUTINE IterativeSolverOption
    END INTERFACE
    CHARACTER(kind = c_char), DIMENSION(:), ALLOCATABLE :: keyc, valc
    ALLOCATE(keyc(LEN(key)+1))
    ALLOCATE(valc(LEN(val)+1))
    CALL c_string_from_f(key, keyc)
    CALL c_string_from_f(val, valc)
    CALL IterativeSolverOption(keyc, valc)
  END SUBROUTINE Iterative_Solver_Option

  !> \brief the lowest eigenvalues of the reduced problem, for the number of roots sought.
  FUNCTION Iterative_Solver_Eigenvalues()
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Iterative_Solver_Eigenvalues
    INTERFACE
      SUBROUTINE IterativeSolverEigenvalues(eigenvalues) BIND(C, name = 'IterativeSolverEigenvalues')
        USE iso_c_binding
        REAL(C_double), DIMENSION(*), INTENT(inout) :: eigenvalues
      END SUBROUTINE IterativeSolverEigenvalues
    END INTERFACE
    ALLOCATE (Iterative_Solver_Eigenvalues(m_nroot))
    CALL IterativeSolverEigenvalues(Iterative_Solver_Eigenvalues)
  END FUNCTION Iterative_Solver_Eigenvalues

  !> @brief Convert from Fortran string to C string
  SUBROUTINE c_string_from_f(fstring, cstring)
    CHARACTER(kind = c_char), DIMENSION(*) :: cstring !< A C char[] big enough to hold the result. No checks are made for overflow.
    CHARACTER(*), INTENT(in) :: fstring !< The fortran string to be converted
    DO i = 1, len_TRIM(fstring)
      cstring(i) = fstring(i : i)
    END DO
    cstring(len_TRIM(fstring) + 1) = C_NULL_CHAR
  END SUBROUTINE c_string_from_f

  !!> Unit testing of IterativeSolver Fortran binding
  SUBROUTINE Iterative_Solver_Test() BIND(C, name = 'IterativeSolverFTest')
    INTEGER, PARAMETER :: n = 100, nroot = 2, nPmax = 20
    DOUBLE PRECISION, DIMENSION (n, n) :: m
    DOUBLE PRECISION, DIMENSION (n, nroot) :: c, g
    DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: p
    DOUBLE PRECISION, DIMENSION (nroot) :: e, error
    INTEGER, DIMENSION(0 : nPmax) :: offsets
    INTEGER, DIMENSION(nPmax) :: indices
    DOUBLE PRECISION, DIMENSION(nPmax) :: coefficients
    DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: pp
    INTEGER :: i, j, root
    DOUBLE PRECISION :: alpha, anharmonicity, threshold
    LOGICAL :: orthogonalize, update
    PRINT *, 'Test Fortran binding of IterativeSolver'
    m = 1
    DO i = 1, n
      m(i, i) = 3 * i
    END DO

    DO irep = 1, 1
      orthogonalize = irep==1
      WRITE (6, *) 'Without P-space, dimension=', n, ', roots=', nroot, ' orthogonalize=', orthogonalize
      threshold = 1d-6; if (orthogonalize) threshold=1d-9
      CALL Iterative_Solver_Linear_Eigensystem_Initialize(n, nroot, &
          thresh = threshold, verbosity = 1, orthogonalize = orthogonalize)
      CALL Iterative_Solver_Option("convergence", "residual")
      c = 0; DO i = 1, nroot; c(i, i) = 1;
      ENDDO
      g = MATMUL(m, c)
      IF(Iterative_Solver_Add_Vector(c, g, p)) THEN
        e = Iterative_Solver_Eigenvalues()
        DO root = 1, nroot
          DO j = 1, n
            c(j, root) = c(j, root) - g(j, root) / (m(j, j) - e(root) + 1d-15)
          END DO
        END DO
      END IF
        IF (Iterative_Solver_End_Iteration(c, g, error)) EXIT
      CALL Iterative_Solver_Finalize
      !write (6,*) 'end of irep loop ',irep
    ENDDO

    DO np = nroot, nPmax, nroot
      ALLOCATE(pp(np, np))
      ALLOCATE(p(np, nroot))
      WRITE (6, *) 'P-space=', nP, ', dimension=', n, ', roots=', nroot
      CALL Iterative_Solver_Linear_Eigensystem_Initialize(n, nroot, thresh = 1d-8, verbosity = 1)
      update = .TRUE.
      offsets(0) = 0
      DO i = 1, nP
        offsets(i) = i
        indices(i) = i
        coefficients(i) = 1
        DO j = 1, nP
          pp(i, j) = m(i, j)
        END DO
      END DO
      c = 0
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
        IF (update) THEN
          DO root = 1, nroot
            DO j = 1, n
              c(j, root) = c(j, root) - g(j, root) / (m(j, j) - e(i) + 1d-15)
            END DO
          END DO
        END IF
        IF (Iterative_Solver_End_Iteration(c, g, error)) EXIT
        g = MATMUL(m, c)
        update = Iterative_Solver_Add_Vector(c, g, p)
      END DO
      CALL Iterative_Solver_Finalize
      DEALLOCATE(p)
      DEALLOCATE(pp)
    END DO

    alpha = 1
    anharmonicity = .5
    WRITE (6, *) 'DIIS, dimension=', n
    CALL Iterative_Solver_DIIS_Initialize(n, thresh = 1d-10, verbosity = 1)
    c = 0;  c(1, 1) = 1
    DO iter = 1, 1000
      DO i = 1, n
        g(i, 1) = (alpha * (i) + anharmonicity * c(i, 1)) * c(i, 1);
        DO j = 1, n
          g(i, 1) = g(i, 1) + (i + j - 2) * c(j, 1);
        END DO
      END DO
      !WRITE (6,*) 'c ',c(:,1)
      !WRITE (6,*) 'g ',g(:,1)
      IF (Iterative_Solver_Add_Vector(c, g, p)) THEN
        DO j = 1, n
          c(j, 1) = c(j, 1) - g(j, 1) / (alpha * (j))
        END DO
      END IF
      IF (Iterative_Solver_End_Iteration(c, g, error)) EXIT
    END DO
    WRITE (6, *) 'error ', error(1), SQRT(dot_PRODUCT(c(:, 1), c(:, 1)))
    CALL Iterative_Solver_Finalize
  END SUBROUTINE Iterative_Solver_Test
END MODULE Iterative_Solver
