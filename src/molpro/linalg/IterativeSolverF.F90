!> @brief IterativeSolver Fortran binding
MODULE Iterative_Solver
    USE, INTRINSIC :: iso_c_binding
    PUBLIC :: Iterative_Solver_Linear_Eigensystem_Initialize, Iterative_Solver_Finalize
    PUBLIC :: Iterative_Solver_Linear_Eigensystem_Initialize_Ranges
    PUBLIC :: Iterative_Solver_DIIS_Initialize, Iterative_Solver_Linear_Equations_Initialize
    PUBLIC :: Iterative_Solver_Optimize_Initialize
    PUBLIC :: Iterative_Solver_Add_Value, Iterative_Solver_End_Iteration
    PUBLIC :: Iterative_Solver_Add_Vector
    PUBLIC :: Iterative_Solver_Solution
    PUBLIC :: Iterative_Solver_Add_P, Iterative_Solver_Suggest_P
    PUBLIC :: Iterative_Solver_Errors
    PUBLIC :: Iterative_Solver_Eigenvalues, Iterative_Solver_Working_Set_Eigenvalues
    PUBLIC :: Iterative_Solver_Print_Statistics
    INTEGER, PUBLIC, PARAMETER :: mpicomm_kind = KIND(c_int64_t)
    PUBLIC :: mpicomm_global, set_mpicomm_compute, mpicomm_self, mpicomm_compute
    PRIVATE
    INTEGER(kind = mpicomm_kind), SAVE :: s_mpicomm_compute=-9999999
    INTEGER(c_size_t) :: m_nq, m_nroot
    INTEGER(c_int), parameter :: hermitian_default=0

    INTERFACE
        SUBROUTINE Iterative_Solver_Print_Statistics() BIND (C, name = 'IterativeSolverPrintStatistics')
        END SUBROUTINE Iterative_Solver_Print_Statistics

        FUNCTION mpicomm_self() BIND(C)
            INTEGER, PARAMETER :: mpicomm_kind = KIND(c_int64_t)
            INTEGER(KIND = mpicomm_kind) :: mpicomm_self
        END FUNCTION mpicomm_self
        FUNCTION mpicomm_global() BIND(C)
            INTEGER, PARAMETER :: mpicomm_kind = KIND(c_int64_t)
            INTEGER(KIND = mpicomm_kind) :: mpicomm_global
        END FUNCTION mpicomm_global
    END INTERFACE

    INTERFACE Iterative_Solver_End_Iteration
        MODULE PROCEDURE Iterative_Solver_End_Iteration1, Iterative_Solver_End_Iteration2
    END INTERFACE Iterative_Solver_End_Iteration

CONTAINS

    FUNCTION mpicomm_compute()
        INTEGER(KIND = mpicomm_kind) :: mpicomm_compute
        if (s_mpicomm_compute .EQ. -9999999) s_mpicomm_compute = mpicomm_global()
        mpicomm_compute = s_mpicomm_compute
    END FUNCTION mpicomm_compute

    SUBROUTINE set_mpicomm_compute(comm)
        INTEGER(KIND = mpicomm_kind), intent(IN) :: comm
        s_mpicomm_compute = comm
    END SUBROUTINE set_mpicomm_compute


    !> \brief Finds the lowest eigensolutions of a matrix. The default algorithm is Davidson's method, i.e. preconditioned Lanczos.
    !> Example of simplest use: @include LinearEigensystemExampleF.F90
    !> Example including use of P space: @include LinearEigensystemExampleF-Pspace-mpi.F90
    SUBROUTINE Iterative_Solver_Linear_Eigensystem_Initialize(nq, nroot, thresh, thresh_value, hermitian, &
        verbosity, pname, mpicomm, algorithm, range)
        INTEGER, INTENT(in) :: nq !< dimension of matrix
        INTEGER, INTENT(in) :: nroot !< number of eigensolutions desired
        DOUBLE PRECISION, INTENT(in), OPTIONAL :: thresh !< Convergence threshold
        DOUBLE PRECISION, INTENT(in), OPTIONAL :: thresh_value !< Value convergence threshold
        LOGICAL, INTENT(in), OPTIONAL :: hermitian !< Whether the underlying kernel is hermitian
        INTEGER, INTENT(in), OPTIONAL :: verbosity !< Print level
        !< One gives a single progress-report line each iteration.
        CHARACTER(len = *), INTENT(in), OPTIONAL :: pname !< Profiler object name
        INTEGER(KIND=mpicomm_kind), INTENT(in), OPTIONAL :: mpicomm !< MPI communicator
        CHARACTER(len = *), INTENT(in), OPTIONAL :: algorithm !< algorithm, eg Davidson
        INTEGER, DIMENSION(2), INTENT(inout), OPTIONAL :: range !< distributed array local range start and end indices
        INTERFACE
            SUBROUTINE Iterative_Solver_Linear_Eigensystem_InitializeC(nq, nroot, range_begin, range_end, thresh, thresh_value, &
                    hermitian, verbosity, pname, mpicomm, algorithm &
                ) BIND(C, name = 'IterativeSolverLinearEigensystemInitialize')
                USE iso_c_binding
                INTEGER(C_size_t), INTENT(in), VALUE :: nq
                INTEGER(C_size_t), INTENT(in), VALUE :: nroot
                INTEGER(C_size_t), INTENT(inout) :: range_begin, range_end
                REAL(c_double), INTENT(in), VALUE :: thresh
                REAL(c_double), INTENT(in), VALUE :: thresh_value
                INTEGER(C_int), INTENT(in), VALUE :: hermitian, verbosity
                CHARACTER(kind = c_char), DIMENSION(*), INTENT(in) :: pname
                INTEGER(C_int64_t), INTENT(in), VALUE :: mpicomm
                CHARACTER(kind = c_char), DIMENSION(*), INTENT(in) :: algorithm
            END SUBROUTINE Iterative_Solver_Linear_Eigensystem_InitializeC
        END INTERFACE
        INTEGER(c_size_t) :: m_range_begin, m_range_end
        INTEGER(c_int) :: hermitianC
        INTEGER(c_int) :: verbosityC
        REAL(c_double) :: threshC
        REAL(c_double) :: thresh_valueC
        CHARACTER(kind = c_char), DIMENSION(:), ALLOCATABLE :: pnameC
        INTEGER(c_int64_t) :: mpicommC
        CHARACTER(kind = c_char), DIMENSION(:), ALLOCATABLE :: algorithmC
        verbosityC = 0
        threshC = 1d-10
        thresh_valueC = 1d50
        hermitianC = hermitian_default
        IF (PRESENT(range)) THEN
            m_range_begin = INT(range(1), kind = c_size_t)
            m_range_end = INT(range(2), kind = c_size_t)
        END IF
        IF (PRESENT(pname)) THEN
            ALLOCATE(pnameC(LEN(pname) + 1))
            CALL c_string_from_f(pname, pnameC)
        ELSE
            ALLOCATE(pnameC(1))
            pnameC(1) = c_null_char
        ENDIF
        IF (PRESENT(algorithm)) THEN
            ALLOCATE(algorithmC(LEN(algorithm) + 1))
            CALL c_string_from_f(algorithm, algorithmC)
        ELSE
            ALLOCATE(algorithmC(1))
            algorithmC(1) = c_null_char
        ENDIF
        m_nq = INT(nq, kind = c_size_t)
        m_nroot = INT(nroot, kind = c_size_t)
        IF (PRESENT(thresh_value)) THEN
            thresh_valueC = thresh_value
        END IF
        IF (PRESENT(thresh)) THEN
            threshC = thresh
        END IF
        IF (PRESENT(hermitian)) THEN
            IF (hermitian) hermitianC = 1
            IF (.NOT. hermitian) hermitianC = 0
        END IF
        IF (PRESENT(verbosity)) THEN
            verbosityC = INT(verbosity, kind = c_int)
        END IF
        IF (PRESENT(mpicomm)) THEN
            mpicommC = INT(mpicomm, kind = c_int64_t)
        ELSE
            mpicommC = mpicomm_compute()
        ENDIF
        CALL Iterative_Solver_Linear_Eigensystem_InitializeC(m_nq, m_nroot, m_range_begin, m_range_end, &
                threshC, thresh_valueC, &
                hermitianC, verbosityC, pnameC, mpicommC, algorithmC)
        IF (PRESENT(range)) THEN
            range(1) = int(m_range_begin)
            range(2) = int(m_range_end)
        END IF
    END SUBROUTINE Iterative_Solver_Linear_Eigensystem_Initialize
    !
    !> \brief Finds the solutions of linear equation systems using a generalisation of Davidson's method, i.e. preconditioned Lanczos
    !> Example of simplest use: @include LinearEquationsExampleF.F90
    !! Example including use of P space: include LinearEquationsExampleF-Pspace.F90
    SUBROUTINE Iterative_Solver_Linear_Equations_Initialize(nq, nroot, rhs, augmented_hessian, thresh, thresh_value, &
        hermitian, verbosity, pname, mpicomm, algorithm, range)
        INTEGER, INTENT(in) :: nq !< dimension of matrix
        INTEGER, INTENT(in) :: nroot !< number of eigensolutions desired
        DOUBLE PRECISION, INTENT(in), DIMENSION(nq, nroot) :: rhs !< the constant right-hand-side of each equation system
        DOUBLE PRECISION, INTENT(in), OPTIONAL :: augmented_hessian
        !< If zero, solve the inhomogeneous equations unmodified. If 1, solve instead
        !< the augmented hessian problem. Other values scale the augmented hessian damping.
        DOUBLE PRECISION, INTENT(in), OPTIONAL :: thresh !< convergence threshold
        DOUBLE PRECISION, INTENT(in), OPTIONAL :: thresh_value !< value convergence threshold
        LOGICAL, INTENT(in), OPTIONAL :: hermitian !< whether the underlying kernel is hermitian
        INTEGER, INTENT(in), OPTIONAL :: verbosity !< how much to print. Default is zero, which prints nothing except errors.
        !< One gives a single progress-report line each iteration.
        CHARACTER(len = *), INTENT(in), OPTIONAL :: pname !< Profiler object name
        INTEGER, INTENT(in), OPTIONAL :: mpicomm !< MPI communicator
        CHARACTER(len = *), INTENT(in), OPTIONAL :: algorithm !< algorithm
        INTEGER, DIMENSION(2), INTENT(inout), OPTIONAL :: range !< distributed array local range start and end indices
        INTERFACE
            SUBROUTINE Iterative_Solver_Linear_Equations_InitializeC(nq, nroot, range_begin, range_end, rhs, &
                    augmented_hessian, thresh, thresh_value, hermitian, verbosity, pname, mpicomm, algorithm) &
                    BIND(C, name = 'IterativeSolverLinearEquationsInitialize')
                USE iso_c_binding
                INTEGER(C_size_t), INTENT(in), VALUE :: nq
                INTEGER(C_size_t), INTENT(in), VALUE :: nroot
                INTEGER(C_size_t), INTENT(inout) :: range_begin, range_end
                REAL(c_double), INTENT(in), DIMENSION(*) :: rhs
                REAL(c_double), INTENT(in), VALUE :: augmented_hessian
                REAL(c_double), INTENT(in), VALUE :: thresh
                REAL(c_double), INTENT(in), VALUE :: thresh_value
                INTEGER(C_int), INTENT(in), VALUE :: hermitian
                INTEGER(C_int), INTENT(in), VALUE :: verbosity
                CHARACTER(kind = c_char), DIMENSION(*), INTENT(in) :: pname
                INTEGER(C_int64_t), INTENT(in), VALUE :: mpicomm
                CHARACTER(kind = c_char), DIMENSION(*), INTENT(in) :: algorithm
            END SUBROUTINE Iterative_Solver_Linear_Equations_InitializeC
        END INTERFACE
        INTEGER(c_size_t) :: m_range_begin, m_range_end
        INTEGER(c_int) :: hermitianC
        INTEGER(c_int) :: verbosityC = 0
        REAL(c_double) :: threshC, thresh_valueC, augmented_hessianC = 0d0
        CHARACTER(kind = c_char), DIMENSION(:), ALLOCATABLE :: pnameC
        INTEGER(c_int64_t) :: mpicommC
        CHARACTER(kind = c_char), DIMENSION(:), ALLOCATABLE :: algorithmC
        threshC = 1d-10
        thresh_valueC = 1d50
        augmented_hessianC = 0d0
        hermitianC = hermitian_default
        IF (PRESENT(range)) THEN
            m_range_begin = INT(range(1), kind = c_size_t)
            m_range_end = INT(range(2), kind = c_size_t)
        END IF
        IF (PRESENT(pname)) THEN
            ALLOCATE(pnameC(LEN(pname) + 1))
            CALL c_string_from_f(pname, pnameC)
        ELSE
            ALLOCATE(pnameC(1))
            pnameC(1) = c_null_char
        ENDIF
        IF (PRESENT(algorithm)) THEN
            ALLOCATE(algorithmC(LEN(algorithm) + 1))
            CALL c_string_from_f(algorithm, algorithmC)
        ELSE
            ALLOCATE(algorithmC(1))
            algorithmC(1) = c_null_char
        ENDIF
        m_nq = INT(nq, kind = c_size_t)
        m_nroot = INT(nroot, kind = c_size_t)
        IF (PRESENT(augmented_hessian)) THEN
            augmented_hessianC = augmented_hessian
        END IF
        IF (PRESENT(thresh)) THEN
            threshC = thresh
        END IF
        IF (PRESENT(thresh_value)) THEN
            thresh_valueC = thresh_value
        END IF
        IF (PRESENT(hermitian)) THEN
            IF (hermitian) hermitianC = 1
            IF (.NOT. hermitian) hermitianC = 0
        END IF
        IF (PRESENT(verbosity)) THEN
            verbosityC = INT(verbosity, kind = c_int)
        END IF
        IF (PRESENT(mpicomm)) THEN
            mpicommC = INT(mpicomm, kind = c_int64_t)
        ELSE
            mpicommC = mpicomm_compute()
        ENDIF
        CALL Iterative_Solver_Linear_Equations_InitializeC(m_nq, m_nroot, m_range_begin, m_range_end, &
                rhs, augmented_hessianC, threshC, thresh_valueC, hermitianC, verbosityC, pnameC, mpicommC, &
            algorithmC)
        IF (PRESENT(range)) THEN
            range(1) = int(m_range_begin)
            range(2) = int(m_range_end)
        END IF
    END SUBROUTINE Iterative_Solver_Linear_Equations_Initialize

    !> \brief Optimization
    !> through the L-BFGS or related methods.
    !> Example of simplest use: @include OptimizeExampleF.F90
    SUBROUTINE Iterative_Solver_Optimize_Initialize(nq, thresh, verbosity, minimize, pname, mpicomm, algorithm, range, thresh_value)
        INTEGER, INTENT(in) :: nq !< dimension of parameter space
        DOUBLE PRECISION, INTENT(in), OPTIONAL :: thresh !< convergence threshold
        INTEGER, INTENT(in), OPTIONAL :: verbosity !< how much to print. Default is zero, which prints nothing except errors.
        !< One gives a single progress-report line each iteration.
        LOGICAL, INTENT(in), OPTIONAL :: minimize !< whether to minimize (default) or maximize
        CHARACTER(len = *), INTENT(in), OPTIONAL :: pname !< Profiler object name
        INTEGER, INTENT(in), OPTIONAL :: mpicomm !< MPI communicator
        CHARACTER(*), INTENT(in), OPTIONAL :: algorithm !< keyword specifying optimization algorithm
        INTEGER, DIMENSION(2), INTENT(inout), OPTIONAL :: range !< distributed array local range start and end indices
        DOUBLE PRECISION, INTENT(in), OPTIONAL :: thresh_value !< convergence threshold for function value
        INTERFACE
!            extern "C" void IterativeSolverOptimizeInitialize(size_t n, size_t* range_begin, size_t* range_end, double thresh,
!                double thresh_value, int verbosity, int minimize, const char* fname,
!                    int64_t fcomm, const char* algorithm);
            SUBROUTINE Iterative_Solver_Optimize_InitializeC(nq, range_begin, range_end, thresh, thresh_value, verbosity, &
                    minimize, pname, mpicomm, algorithm) BIND(C, name = 'IterativeSolverOptimizeInitialize')
                USE iso_c_binding
                INTEGER(C_size_t), INTENT(in), VALUE :: nq
                INTEGER(C_size_t), INTENT(inout) :: range_begin, range_end
                REAL(c_double), INTENT(in), VALUE :: thresh, thresh_value
                INTEGER(C_int), INTENT(in), VALUE :: verbosity
                INTEGER(C_int), INTENT(in), VALUE :: minimize
                CHARACTER(kind = c_char), DIMENSION(*), INTENT(in) :: algorithm
                CHARACTER(kind = c_char), DIMENSION(*), INTENT(in) :: pname
                INTEGER(C_int64_t), INTENT(in), VALUE :: mpicomm
            END SUBROUTINE Iterative_Solver_Optimize_InitializeC
        END INTERFACE
        INTEGER(c_size_t) :: m_range_begin, m_range_end
        INTEGER(c_int) :: verbosityC, minimizeC
        REAL(c_double) :: threshC, thresh_valueC
        CHARACTER(kind = c_char), DIMENSION(:), ALLOCATABLE :: pnameC
        INTEGER(c_int64_t) :: mpicommC
        CHARACTER(kind = c_char), DIMENSION(:), ALLOCATABLE :: algorithmC
        IF (PRESENT(range)) THEN
            m_range_begin = INT(range(1), kind = c_size_t)
            m_range_end = INT(range(2), kind = c_size_t)
        END IF
        IF (PRESENT(pname)) THEN
            ALLOCATE(pnameC(LEN(pname) + 1))
            CALL c_string_from_f(pname, pnameC)
        ELSE
            ALLOCATE(pnameC(1))
            pnameC(1) = c_null_char
        ENDIF
        m_nq = INT(nq, kind = c_size_t)
        m_nroot = 1
        IF (PRESENT(thresh)) THEN
            threshC = thresh
        ELSE
            threshC = 1d-10
        END IF
        IF (PRESENT(thresh_value)) THEN
            thresh_valueC = thresh_value
        ELSE
            thresh_valueC = 1d50
        END IF
        IF (PRESENT(verbosity)) THEN
            verbosityC = INT(verbosity, kind = c_int)
        ELSE
            verbosityC = 0
        END IF
        IF (PRESENT(algorithm)) THEN
            ALLOCATE(algorithmC(LEN(algorithm) + 1))
            CALL c_string_from_f(algorithm, algorithmC)
        ELSE
            ALLOCATE(algorithmC(1))
            algorithmC(1) = c_null_char
        ENDIF
        minimize_C = 1
        IF (PRESENT(minimize)) THEN
            IF (.NOT. minimize) minimizeC = 0
        END IF
        IF (PRESENT(mpicomm)) THEN
            mpicommC = INT(mpicomm, kind = c_int64_t)
        ELSE
            mpicommC = mpicomm_compute()
        ENDIF
        CALL Iterative_Solver_Optimize_InitializeC(m_nq, m_range_begin, m_range_end, threshC, thresh_valueC, verbosityC, &
                minimizeC, pnameC, mpicommC, algorithmC)
        IF (PRESENT(range)) THEN
            range(1) = int(m_range_begin)
            range(2) = int(m_range_end)
        END IF
    END SUBROUTINE Iterative_Solver_Optimize_Initialize

    !> \brief Accelerated convergence of non-linear equations
    !> through the DIIS or related methods.
    !> Example of simplest use: @include DIISExampleF.F90
    SUBROUTINE Iterative_Solver_DIIS_Initialize(nq, thresh, verbosity, pname, mpicomm, algorithm, range)
        INTEGER, INTENT(in) :: nq !< dimension of parameter space
        DOUBLE PRECISION, INTENT(in), OPTIONAL :: thresh !< convergence threshold
        INTEGER, INTENT(in), OPTIONAL :: verbosity !< how much to print. Default is zero, which prints nothing except errors.
        !< One gives a single progress-report line each iteration.
        CHARACTER(len = *), INTENT(in), OPTIONAL :: pname !< Profiler object name
        INTEGER, INTENT(in), OPTIONAL :: mpicomm !< Profiler communicator
        CHARACTER(len = *), INTENT(in), OPTIONAL :: algorithm !< algorithm, eg DIIS
        INTEGER, DIMENSION(2), INTENT(inout), OPTIONAL :: range !< distributed array local range start and end indices
        INTERFACE
            SUBROUTINE Iterative_Solver_DIIS_InitializeC(nq, range_begin, range_end, thresh, verbosity, &
                    pname, mpicomm, algorithm) BIND(C, name = 'IterativeSolverNonLinearEquationsInitialize')
                USE iso_c_binding
                INTEGER(C_size_t), INTENT(in), VALUE :: nq
                INTEGER(C_size_t), INTENT(inout) :: range_begin, range_end
                REAL(c_double), INTENT(in), VALUE :: thresh
                INTEGER(C_int), INTENT(in), VALUE :: verbosity
                CHARACTER(kind = c_char), DIMENSION(*), INTENT(in) :: pname
                INTEGER(C_int64_t), INTENT(in), VALUE :: mpicomm
                CHARACTER(kind = c_char), DIMENSION(*), INTENT(in) :: algorithm
            END SUBROUTINE Iterative_Solver_DIIS_InitializeC
        END INTERFACE
        INTEGER(c_size_t) :: m_range_begin, m_range_end
        INTEGER(c_int) :: verbosityC
        REAL(c_double) :: threshC
        CHARACTER(kind = c_char), DIMENSION(:), ALLOCATABLE :: pnameC
        INTEGER(c_int64_t) :: mpicommC
        CHARACTER(kind = c_char), DIMENSION(:), ALLOCATABLE :: algorithmC
        IF (PRESENT(range)) THEN
            m_range_begin = INT(range(1), kind = c_size_t)
            m_range_end = INT(range(2), kind = c_size_t)
        END IF
        IF (PRESENT(pname)) THEN
            ALLOCATE(pnameC(LEN(pname) + 1))
            CALL c_string_from_f(pname, pnameC)
        ELSE
            ALLOCATE(pnameC(1))
            pnameC(1) = c_null_char
        ENDIF
        IF (PRESENT(algorithm)) THEN
            ALLOCATE(algorithmC(LEN(algorithm) + 1))
            CALL c_string_from_f(algorithm, algorithmC)
        ELSE
            ALLOCATE(algorithmC(1))
            algorithmC(1) = c_null_char
        ENDIF
        m_nq = INT(nq, kind = c_size_t)
        m_nroot = 1
        IF (PRESENT(thresh)) THEN
            threshC = thresh
        ELSE
            threshC = 1d-10
        END IF
        IF (PRESENT(verbosity)) THEN
            verbosityC = INT(verbosity, kind = c_int)
        ELSE
            verbosityC = 0
        END IF
        IF (PRESENT(mpicomm)) THEN
            mpicommC = INT(mpicomm, kind = c_int64_t)
        ELSE
            mpicommC = mpicomm_compute()
        ENDIF
        CALL Iterative_Solver_DIIS_InitializeC(m_nq, m_range_begin, m_range_end, threshC, verbosityC, &
                pnameC, mpicommC, algorithm)
        IF (PRESENT(range)) THEN
            range(1) = int(m_range_begin)
            range(2) = int(m_range_end)
        END IF
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

    !> \brief Take a current solution, value and residual, add it to the expansion set, and return working-set residual.
    !> \param value On input, the current objective function value.
    !> \param parameters On input, the current solution or expansion vector. On exit, undefined.
    !> \param action On input, the residual for parameters.
    !> On exit, the expected residual of the new working set.
    !> \param synchronize Whether to synchronize any distributed storage of parameters and action before return.
    !>        Unnecessary if the client preconditioner is diagonal, but otherwise should be done.
    !>        The default is the safe .TRUE. but can be .FALSE. if appropriate.
    !> \return whether it is expected that the client should make an update, based on the returned parameters and residual, before
    !> the subsequent call to Iterative_Solver_End_Iteration()
    FUNCTION Iterative_Solver_Add_Value(value, parameters, action, synchronize)
        USE iso_c_binding
        LOGICAL :: Iterative_Solver_Add_Value
        DOUBLE PRECISION, INTENT(in) :: value
        DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: parameters
        DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: action
        LOGICAL, INTENT(in), OPTIONAL :: synchronize
        INTERFACE
            FUNCTION Iterative_Solver_Add_Value_C(value, parameters, action, lsync) &
                    BIND(C, name = 'IterativeSolverAddValue')
                USE iso_c_binding
                INTEGER(c_size_t) Iterative_Solver_Add_Value_C
                REAL(c_double), VALUE, INTENT(in) :: value
                REAL(c_double), DIMENSION(*), INTENT(inout) :: parameters
                REAL(c_double), DIMENSION(*), INTENT(inout) :: action
                INTEGER(c_int), INTENT(in), VALUE :: lsync
            END FUNCTION Iterative_Solver_Add_Value_C
        END INTERFACE
        INTEGER(c_int) :: lsyncC
        lsyncC = 1
        IF (PRESENT(synchronize)) THEN
            IF (.NOT. synchronize) lsyncC = 0
        END IF
        Iterative_Solver_Add_Value = Iterative_Solver_Add_Value_C(value, parameters, action, lsyncC).NE.0
    END FUNCTION Iterative_Solver_Add_Value
    !
    !
    !> \brief Take, typically, a current solution and residual, add it to the expansion set, and return new solution.
    !> In the context of Lanczos-like linear methods, the input will be a current expansion vector and the result of
    !> acting on it with the matrix, and the output will be a new expansion vector.
    !> \param parameters On input, the current solution or expansion vector. On exit, the interpolated solution vector.
    !> \param action On input, the residual for parameters (non-linear), or action of matrix on parameters (linear).
    !> On exit, the expected (non-linear) or actual (linear) residual of the interpolated parameters.
    !> \param synchronize Whether to synchronize any distributed storage of parameters and action before return. Unnecessary if the
    !> client preconditioner is diagonal, but otherwise should be done. The default is the safe .TRUE. but can be .FALSE. if
    !> appropriate.
    !> \return whether it is expected that the client should make an update, based on the returned parameters and residual, before
    !> the subsequent call to Iterative_Solver_End_Iteration()

    FUNCTION Iterative_Solver_Add_Vector(parameters, action, synchronize)
        USE iso_c_binding
        INTEGER :: Iterative_Solver_Add_Vector
        DOUBLE PRECISION, DIMENSION(:,:), INTENT(inout) :: parameters
        DOUBLE PRECISION, DIMENSION(:,:), INTENT(inout) :: action
        LOGICAL, INTENT(in), OPTIONAL :: synchronize
        INTERFACE
            FUNCTION Add_Vector_C(buffer_size, parameters, action, lsync) &
                    BIND(C, name = 'IterativeSolverAddVector')
                USE, INTRINSIC :: iso_c_binding
                INTEGER(c_size_t) Add_Vector_C
                INTEGER(c_size_t), INTENT(in), VALUE :: buffer_size
                REAL(c_double), DIMENSION(*), INTENT(inout) :: parameters
                REAL(c_double), DIMENSION(*), INTENT(inout) :: action
                INTEGER(c_int), INTENT(in), VALUE :: lsync
            END FUNCTION Add_Vector_C
        END INTERFACE
        DOUBLE PRECISION, DIMENSION(0) :: pdummy
        INTEGER(c_int) :: lsyncC
        lsyncC = 1
        IF (PRESENT(synchronize)) THEN
            IF (.NOT. synchronize) lsyncC = 0
        END IF
        Iterative_Solver_Add_Vector = int(&
            Add_Vector_C(int(ubound(parameters, 2) - lbound(parameters, 2) + 1, c_size_t), &
                parameters, action, lsyncC))
    END FUNCTION Iterative_Solver_Add_Vector
    !
    SUBROUTINE Iterative_Solver_Solution(roots, parameters, action, synchronize)
        USE iso_c_binding
        INTEGER, INTENT(in), DIMENSION(:) :: roots  !< Array containing root indices
        DOUBLE PRECISION, DIMENSION(:,:), INTENT(inout) :: parameters
        DOUBLE PRECISION, DIMENSION(:,:), INTENT(inout) :: action
        LOGICAL, INTENT(in), OPTIONAL :: synchronize
        INTERFACE
            SUBROUTINE Solution_C(nroot, roots, parameters, action, lsync) &
                    BIND(C, name = 'IterativeSolverSolution')
                USE iso_c_binding
                INTEGER(c_int), VALUE :: nroot
                INTEGER(c_int), INTENT(in), DIMENSION(nroot) :: roots
                REAL(c_double), DIMENSION(*), INTENT(inout) :: parameters
                REAL(c_double), DIMENSION(*), INTENT(inout) :: action
                INTEGER(c_int), INTENT(in), VALUE :: lsync
            END SUBROUTINE
        END INTERFACE
        INTEGER(c_int), DIMENSION(SIZE(roots)) :: rootsC
        INTEGER(c_int) :: nroot
        INTEGER(c_int) :: lsyncC
        lsyncC = 1
        IF (PRESENT(synchronize)) THEN
            IF (.NOT. synchronize) lsyncC = 0
        END IF
        nroot = INT(size(roots), c_int)
        DO i = 1, size(roots)
            rootsC(i) = INT(roots(i)-1, kind = c_int)
        ENDDO
        call Solution_C(nroot, rootsC, parameters, action, lsyncC)
    END SUBROUTINE
    !
    !>@brief For most methods, does nothing; for Optimize it is required.
    !> Also write progress to standard output
    !> \param solution On exit, the new vectors targeting roots in the working set.
    !> \param residual On entry, the preconditioned residual. On exit, undefined.
    !> \param synchronize Whether to synchronize any distributed storage of parameters before return.
    !>        The default is the safe .TRUE. but can be .FALSE. if appropriate.
    !> \return The size of the working set
    FUNCTION Iterative_Solver_End_Iteration1(solution, residual, synchronize, buffer_size)
        USE iso_c_binding
        INTEGER :: Iterative_Solver_End_Iteration1
        DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: solution
        DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: residual
        LOGICAL, INTENT(in), OPTIONAL :: synchronize
        INTEGER, INTENT(in), OPTIONAL :: buffer_size
    INTERFACE
            FUNCTION Iterative_Solver_End_Iteration_C(buffer_size, solution, residual, lsync) &
                BIND(C, name = 'IterativeSolverEndIteration')
                USE iso_c_binding
                INTEGER(c_size_t) Iterative_Solver_End_Iteration_C
                INTEGER(c_size_t), INTENT(in), VALUE :: buffer_size
                REAL(c_double), DIMENSION(*), INTENT(inout) :: solution
                REAL(c_double), DIMENSION(*), INTENT(inout) :: residual
                INTEGER(c_int), INTENT(in), VALUE :: lsync
            END FUNCTION Iterative_Solver_End_Iteration_C
        END INTERFACE
        INTEGER(c_int) :: lsyncC
        INTEGER(c_size_t) :: buffer_sizeC
        buffer_sizeC = 1
        IF (PRESENT(buffer_size)) buffer_sizeC = buffer_size
        lsyncC = 1
        IF (PRESENT(synchronize)) THEN
            IF (.NOT. synchronize) lsyncC = 0
        END IF
        Iterative_Solver_End_Iteration1 = int(Iterative_Solver_End_Iteration_C( &
            buffer_sizeC, &
            solution, residual, lsyncC))
    END FUNCTION Iterative_Solver_End_Iteration1

    FUNCTION Iterative_Solver_End_Iteration2(solution, residual, synchronize)
        USE iso_c_binding
        INTEGER :: Iterative_Solver_End_Iteration2
        DOUBLE PRECISION, DIMENSION(:,:), INTENT(inout) :: solution
        DOUBLE PRECISION, DIMENSION(:,:), INTENT(inout) :: residual
        LOGICAL, INTENT(in), OPTIONAL :: synchronize
        Iterative_Solver_End_Iteration2 = Iterative_Solver_End_Iteration1( &
            solution, residual, synchronize, &
        buffer_size = ubound(solution, 2) - lbound(solution, 2) + 1 &
        )
    END FUNCTION Iterative_Solver_End_Iteration2


    !> \brief add P-space vectors to the expansion set, and return new solution.
    !> \param nP the number of P-space vectors to add
    !> \param offsets specifies the start point in indices and coefficients that defines each vector.
    !> \param indices Index in the full space of a contribution to a new P vector
    !> \param coefficients Value of a contribution to a new P vector
    !> \param pp The P-P block of the matrix, dimensioned (number of existing P + nP, nP)
    !> \param parameters On input, the current solution or expansion vector. On exit, the interpolated solution vector.
    !> \param action On input, the residual for parameters (non-linear), or action of matrix on parameters (linear).
    !> \param fproc
    !> \param synchronize Whether to synchronize any distributed storage of parameters and action before return.
    !>        Unnecessary if the client preconditioner is diagonal, but otherwise should be done.
    !>        The default is the safe .TRUE. but can be .FALSE. if appropriate.
    !> On exit, the expected (non-linear) or actual (linear) residual of the interpolated parameters.
    FUNCTION Iterative_Solver_Add_P(nP, offsets, indices, coefficients, pp, parameters, action, fproc, synchronize)
        USE iso_c_binding
        INTEGER :: Iterative_Solver_Add_P
        INTEGER, INTENT(in) :: nP
        INTEGER, INTENT(in), DIMENSION(0:nP) :: offsets
        INTEGER, INTENT(in), DIMENSION(offsets(nP)) :: indices
        DOUBLE PRECISION, DIMENSION(offsets(nP)), INTENT(in) :: coefficients
        DOUBLE PRECISION, DIMENSION(*), INTENT(in) :: pp
        DOUBLE PRECISION, DIMENSION(:,:), INTENT(inout) :: parameters
        DOUBLE PRECISION, DIMENSION(:,:), INTENT(inout) :: action
        LOGICAL, INTENT(in), OPTIONAL :: synchronize
        EXTERNAL fproc
        INTERFACE
            FUNCTION IterativeSolverAddPC(buffer_size, nP, offsets, indices, coefficients, pp, parameters, action, &
                    lsync, func) BIND(C, name = 'IterativeSolverAddP')
                USE, INTRINSIC :: iso_c_binding
                INTEGER(c_size_t) IterativeSolverAddPC
                INTEGER(c_size_t), INTENT(in), VALUE :: buffer_size
                INTEGER(c_size_t), INTENT(in), VALUE :: nP
                INTEGER(c_size_t), INTENT(in), DIMENSION(0:nP) :: offsets
                INTEGER(c_size_t), INTENT(in), DIMENSION(offsets(nP)) :: indices
                REAL(c_double), DIMENSION(offsets(nP)), INTENT(in) :: coefficients
                REAL(c_double), DIMENSION(*), INTENT(in) :: pp
                REAL(c_double), DIMENSION(*), INTENT(inout) :: parameters
                REAL(c_double), DIMENSION(*), INTENT(inout) :: action
                INTEGER(c_int), INTENT(in), VALUE :: lsync
                TYPE(C_FUNPTR), INTENT(IN), VALUE :: func
            END FUNCTION IterativeSolverAddPC
        END INTERFACE
        INTEGER(c_size_t), DIMENSION(0:nP) :: offsetsC
        INTEGER(c_size_t), DIMENSION(SIZE(indices)) :: indicesC
        INTEGER(c_int) :: lsyncC
        TYPE(C_FUNPTR) :: cproc
        cproc = C_FUNLOC(fproc)
        lsyncC = 1
        IF (PRESENT(synchronize)) THEN
            IF (.NOT. synchronize) lsyncC = 0
        END IF
        offsetsC = INT(offsets, c_size_t)
        do i = 1, offsets(nP)
            indicesC(i) = INT(indices(i) - 1, c_size_t) ! 1-base to 0-base
        end do
        Iterative_Solver_Add_P = int(IterativeSolverAddPC( &
            INT(ubound(parameters,2)-lbound(parameters,2)+1, c_size_t), &
            INT(nP, c_size_t), offsetsC, indicesC, coefficients, &
                pp, parameters, action, lsyncC, cproc))
    END FUNCTION Iterative_Solver_Add_P

    !> \brief Take an existing solution and its residual, and suggest P vectors
    !> \param solution On input, the current solution.
    !> \param residual On input, the residual for solution.
    !> \param indices On exit, the most important base vectors
    !> \param threshold IterativeSolver vectors whose predicted contribution is less than
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
        indicesC = 0
        maximumNumber = INT(size(indices), c_size_t)
        !write (6,*) 'fortran suggestP, maximumNumber=',size(indices)
        !write (6,*) 'fortran suggestP, indices=',indices
        IF (PRESENT(threshold)) thresholdC = threshold
        Iterative_Solver_Suggest_P = INT(&
                IterativeSolverSuggestP(solution, residual, maximumNumber, thresholdC, indicesC) &
                )
        do i = 1, Iterative_Solver_Suggest_P
            indices(i) = int(indicesC(i)) + 1
        end do
    END FUNCTION Iterative_Solver_Suggest_P

    !> \brief errors for each root
    FUNCTION Iterative_Solver_Errors()
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Iterative_Solver_Errors
        INTERFACE
            SUBROUTINE IterativeSolverErrors(errors) BIND(C, name = 'IterativeSolverErrors')
                USE iso_c_binding
                REAL(C_double), DIMENSION(*), INTENT(inout) :: errors
            END SUBROUTINE IterativeSolverErrors
        END INTERFACE
        ALLOCATE (Iterative_Solver_Errors(m_nroot))
        CALL IterativeSolverErrors(Iterative_Solver_Errors)
    END FUNCTION Iterative_Solver_Errors

    !> \brief the lowest eigenvalues of the reduced problem, for the number of roots sought.
    FUNCTION Iterative_Solver_Eigenvalues()
        !DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Iterative_Solver_Eigenvalues
        DOUBLE PRECISION, DIMENSION(m_nroot) :: Iterative_Solver_Eigenvalues
        INTERFACE
            SUBROUTINE IterativeSolverEigenvalues(eigenvalues) BIND(C, name = 'IterativeSolverEigenvalues')
                USE iso_c_binding
                REAL(C_double), DIMENSION(*), INTENT(inout) :: eigenvalues
            END SUBROUTINE IterativeSolverEigenvalues
        END INTERFACE
        !ALLOCATE (Iterative_Solver_Eigenvalues(m_nroot))
        CALL IterativeSolverEigenvalues(Iterative_Solver_Eigenvalues)
    END FUNCTION Iterative_Solver_Eigenvalues

    !> \brief the eigenvalues of the reduced problem, for the number of roots in working set (not yet converged).
    FUNCTION Iterative_Solver_Working_Set_Eigenvalues(working_set_size)
        !DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Iterative_Solver_Working_Set_Eigenvalues
        INTEGER, INTENT(in) :: working_set_size
        DOUBLE PRECISION, DIMENSION(working_set_size) :: Iterative_Solver_Working_Set_Eigenvalues
        INTERFACE
            SUBROUTINE IterativeSolverWorkingSetEigenvalues(eigenvalues) BIND(C, name = 'IterativeSolverWorkingSetEigenvalues')
                USE iso_c_binding
                REAL(C_double), DIMENSION(*), INTENT(inout) :: eigenvalues
            END SUBROUTINE IterativeSolverWorkingSetEigenvalues
        END INTERFACE
        !ALLOCATE (Iterative_Solver_Working_Set_Eigenvalues(int(working_set_size, c_size_t)))
        CALL IterativeSolverWorkingSetEigenvalues(Iterative_Solver_Working_Set_Eigenvalues)
    END FUNCTION Iterative_Solver_Working_Set_Eigenvalues
    !> @brief Convert from Fortran string to C string
    SUBROUTINE c_string_from_f(fstring, cstring)
        CHARACTER(kind = c_char), DIMENSION(*) :: cstring !< A C char[] big enough to hold the result. No checks are made for overflow.
        CHARACTER(*), INTENT(in) :: fstring !< The fortran string to be converted
        DO i = 1, len_TRIM(fstring)
            cstring(i) = fstring(i:i)
        END DO
        cstring(len_TRIM(fstring) + 1) = C_NULL_CHAR
    END SUBROUTINE c_string_from_f

    !!> Unit testing of IterativeSolver Fortran binding
    !SUBROUTINE Iterative_Solver_Test() BIND(C, name = 'IterativeSolverFTest')
    !  INTEGER, PARAMETER :: n = 100, nroot = 2, nPmax = 20
    !  DOUBLE PRECISION, DIMENSION (n, n) :: m
    !  DOUBLE PRECISION, DIMENSION (n, nroot) :: c, g
    !  DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: p
    !  DOUBLE PRECISION, DIMENSION (nroot) :: e, error
    !  INTEGER, DIMENSION(0 : nPmax) :: offsets
    !  INTEGER, DIMENSION(nPmax) :: indices
    !  DOUBLE PRECISION, DIMENSION(nPmax) :: coefficients
    !  DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: pp
    !  INTEGER :: i, j, root
    !  DOUBLE PRECISION :: alpha, anharmonicity, threshold
    !  INTEGER ::  working_set_size
    !  PRINT *, 'Test Fortran binding of IterativeSolver'
    !  m = 1
    !  DO i = 1, n
    !    m(i, i) = 3 * i
    !  END DO

    !  DO irep = 1, 1
    !    WRITE (6, *) 'Without P-space, dimension=', n, ', roots=', nroot
    !    threshold = 1d-6;
    !    CALL Iterative_Solver_Linear_Eigensystem_Initialize(n, nroot, &
    !        thresh = threshold, verbosity = 1)
    !    CALL Iterative_Solver_Option("convergence", "residual")
    !    c = 0; DO i = 1, nroot; c(i, i) = 1;
    !    ENDDO
    !    g = MATMUL(m, c)
    !    IF(Add_Vector(c, g, p) .GT. 0) THEN
    !      e = Iterative_Solver_Eigenvalues()
    !      DO root = 1, nroot
    !        DO j = 1, n
    !          c(j, root) = c(j, root) - g(j, root) / (m(j, j) - e(root) + 1d-15)
    !        END DO
    !      END DO
    !    END IF
    !      IF (Iterative_Solver_End_Iteration(c, g, error)) EXIT
    !    CALL Iterative_Solver_Finalize
    !    !write (6,*) 'end of irep loop ',irep
    !  ENDDO

    !  DO np = nroot, nPmax, nroot
    !    ALLOCATE(pp(np, np))
    !    ALLOCATE(p(np, nroot))
    !    WRITE (6, *) 'P-space=', nP, ', dimension=', n, ', roots=', nroot
    !    CALL Iterative_Solver_Linear_Eigensystem_Initialize(n, nroot, thresh = 1d-8, verbosity = 1)
    !    update = nroot
    !    offsets(0) = 0
    !    DO i = 1, nP
    !      offsets(i) = i
    !      indices(i) = i
    !      coefficients(i) = 1
    !      DO j = 1, nP
    !        pp(i, j) = m(i, j)
    !      END DO
    !    END DO
    !    c = 0
    !    CALL Iterative_Solver_Add_P(nP, offsets, indices, coefficients, pp, c, g, p)
    !    DO iter = 1, 10
    !      e = Iterative_Solver_Eigenvalues()
    !      DO root = 1, nroot
    !        DO i = 1, nP
    !          DO j = 1, n
    !            g(j, root) = g(j, root) + m(j, indices(i)) * p(i, root)
    !          END DO
    !        END DO
    !      END DO
    !      IF (update > 0) THEN
    !        DO root = 1, nroot
    !          DO j = 1, n
    !            c(j, root) = c(j, root) - g(j, root) / (m(j, j) - e(i) + 1d-15)
    !          END DO
    !        END DO
    !      END IF
    !      IF (Iterative_Solver_End_Iteration(c, g, error)) EXIT
    !      g = MATMUL(m, c)
    !      update = Add_Vector(c, g, p)
    !    END DO
    !    CALL Iterative_Solver_Finalize
    !    DEALLOCATE(p)
    !    DEALLOCATE(pp)
    !  END DO

    !  stop
    !  alpha = 1
    !  anharmonicity = .5
    !  WRITE (6, *) 'DIIS, dimension=', n
    !  CALL Iterative_Solver_DIIS_Initialize(n, thresh = 1d-10, verbosity = 1)
    !  c = 0;  c(1, 1) = 1
    !  DO iter = 1, 1000
    !    DO i = 1, n
    !      g(i, 1) = (alpha * (i) + anharmonicity * c(i, 1)) * c(i, 1);
    !      DO j = 1, n
    !        g(i, 1) = g(i, 1) + (i + j - 2) * c(j, 1);
    !      END DO
    !    END DO
    !    !WRITE (6,*) 'c ',c(:,1)
    !    !WRITE (6,*) 'g ',g(:,1)
    !    IF (Add_Vector(c, g, p) .GT. 0) THEN
    !      DO j = 1, n
    !        c(j, 1) = c(j, 1) - g(j, 1) / (alpha * (j))
    !      END DO
    !    END IF
    !    IF (Iterative_Solver_End_Iteration(c, g, error)) EXIT
    !  END DO
    !  WRITE (6, *) 'error ', error(1), SQRT(dot_PRODUCT(c(:, 1), c(:, 1)))
    !  CALL Iterative_Solver_Finalize
    !END SUBROUTINE Iterative_Solver_Test
END MODULE Iterative_Solver

!> @examples LinearEigensystemExampleF.F90
!> This is an example of simplest use of the LinearEigensystem framework for iterative
!> finding of the lowest few eigensolutions of a large matrix.

!> @example LinearEquationsExampleF.F90
!> This is an examples of simplest use of the LinearEquations framework for iterative
!> solution of linear equations