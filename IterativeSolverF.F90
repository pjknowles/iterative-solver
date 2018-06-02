!> @brief IterativeSolver Fortran binding
MODULE Iterative_Solver
 USE iso_c_binding
 PUBLIC :: Iterative_Solver_Linear_Eigensystem_Initialize, Iterative_Solver_Finalize
 PUBLIC :: Iterative_Solver_Add_Vector, Iterative_Solver_End_Iteration
 PUBLIC :: Iterative_Solver_Add_P
 PUBLIC :: Iterative_Solver_Eigenvalues
 PRIVATE
 INTEGER(c_size_t), PRIVATE :: m_nq, m_nroot
 
CONTAINS
 
!> \brief Finds the lowest eigensolutions of a matrix using Davidson's method, i.e. preconditioned Lanczos
!> Example of simplest use: @include LinearEigensystemExampleF.F90
 SUBROUTINE Iterative_Solver_Linear_Eigensystem_Initialize(nq,nroot,thresh,maxIterations,verbosity)
  INTEGER, INTENT(in) :: nq !< dimension of matrix
  INTEGER, INTENT(in) :: nroot !< number of eigensolutions desired
  DOUBLE PRECISION, INTENT(in), OPTIONAL :: thresh !< convergence threshold
  INTEGER, INTENT(in), OPTIONAL :: maxIterations !< maximum number of iterations
  INTEGER, INTENT(in), OPTIONAL :: verbosity !< how much to print. Default is zero, which prints nothing except errors
  INTERFACE
   SUBROUTINE Iterative_Solver_Linear_Eigensystem_InitializeC(nq,nroot,thresh,maxIterations,verbosity) &
        BIND(C,name='IterativeSolverLinearEigensystemInitialize')
    USE iso_c_binding
    INTEGER(C_size_t), INTENT(in), VALUE :: nq
    INTEGER(C_size_t), INTENT(in), VALUE :: nroot
    REAL(c_double), INTENT(in), VALUE :: thresh
    INTEGER(C_int), INTENT(in), VALUE :: maxIterations
    INTEGER(C_int), INTENT(in), VALUE :: verbosity
   END SUBROUTINE Iterative_Solver_Linear_Eigensystem_InitializeC
  END INTERFACE
  INTEGER(c_int) :: verbosityC, maxIterationsC
  REAL(c_double) :: threshC
  m_nq=INT(nq,kind=c_size_t)
  m_nroot=INT(nroot,kind=c_size_t)
  IF (PRESENT(thresh)) THEN
   threshC=thresh
  ELSE
   threshC=0
  END IF
  IF (PRESENT(maxIterations)) THEN
   maxIterationsC=INT(maxIterations,kind=c_int)
  ELSE
   maxIterationsC=0
  END IF
  IF (PRESENT(verbosity)) THEN
   verbosityC=INT(verbosity,kind=c_int)
  ELSE
   verbosityC=0
  END IF
  CALL Iterative_Solver_Linear_Eigensystem_InitializeC(m_nq,m_nroot,threshC, maxIterationsC, verbosityC)
 END SUBROUTINE Iterative_Solver_Linear_Eigensystem_Initialize

!> \brief Terminate the iterative solver
 SUBROUTINE Iterative_Solver_Finalize
  INTERFACE
   SUBROUTINE IterativeSolverFinalize() BIND(C,name='IterativeSolverFinalize')
    USE iso_c_binding
   END SUBROUTINE IterativeSolverFinalize
  END INTERFACE
  CALL IterativeSolverFinalize
 END SUBROUTINE Iterative_Solver_Finalize
 
!> \brief Take, typically, a current solution and residual, add it to the expansion set, and return new solution.
!> In the context of Lanczos-like linear methods, the input will be a current expansion vector and the result of
!> acting on it with the matrix, and the output will be a new expansion vector.
!> \param parameters On input, the current solution or expansion vector. On exit, the interpolated solution vector.
!> \param action On input, the residual for parameters (non-linear), or action of matrix on parameters (linear). On exit, the expected (non-linear) or actual (linear) residual of the interpolated parameters.
!> \param parametersP On exit, the interpolated solution projected onto the P space.
 SUBROUTINE Iterative_Solver_Add_Vector(parameters,action,parametersP)
  USE iso_c_binding
  DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: parameters
  DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: action
  DOUBLE PRECISION, DIMENSION(*), INTENT(inout), OPTIONAL :: parametersP
  INTERFACE
   SUBROUTINE Iterative_Solver_Add_Vector_C(parameters,action,parametersP) &
        BIND(C,name='IterativeSolverAddVector')
    USE iso_c_binding
    REAL(c_double), DIMENSION(*), INTENT(inout) :: parameters
    REAL(c_double), DIMENSION(*), INTENT(inout) :: action
    REAL(c_double), DIMENSION(*), INTENT(inout) :: parametersP
   END SUBROUTINE Iterative_Solver_Add_Vector_C
  END INTERFACE
  DOUBLE PRECISION, DIMENSION(0) :: pdummy
  IF (PRESENT(parametersP)) THEN
   CALL Iterative_Solver_Add_Vector_C(parameters,action,parametersP)
  ELSE
   CALL Iterative_Solver_Add_Vector_C(parameters,action,pdummy)
  END IF
 END SUBROUTINE Iterative_Solver_Add_Vector
 
!>@brief Take the updated solution vector set, and adjust it if necessary so that it becomes the vector to
!> be used in the next iteration; this is done only in the case of linear solvers where the orthogonalize option is set.
!> Also calculate the degree of convergence, and write progress to standard output
!> \param solution The current solution, after interpolation and updating with the preconditioned residual.
!> \param residual The residual after interpolation.
!> \param error Error indicator for each sought root.
!> \return .TRUE. if convergence reached for all roots
 LOGICAL FUNCTION Iterative_Solver_End_Iteration(solution,residual,error)
  USE iso_c_binding
  DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: solution
  DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: residual
  DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: error
  INTERFACE
   INTEGER(c_int) FUNCTION Iterative_Solver_End_Iteration_C(solution,residual,error) &
        BIND(C,name='IterativeSolverEndIteration')
    USE iso_c_binding
    REAL(c_double), DIMENSION(*), INTENT(inout) :: solution
    REAL(c_double), DIMENSION(*), INTENT(inout) :: residual
    REAL(c_double), DIMENSION(*), INTENT(inout) :: error
   END FUNCTION Iterative_Solver_End_Iteration_C
  END INTERFACE
  Iterative_Solver_End_Iteration = &
       Iterative_Solver_End_Iteration_C(solution,residual,error).NE.0
 END FUNCTION Iterative_Solver_End_Iteration
 
 
!> \brief add P-space vectors to the expansion set, and return new solution.
!> \param nP the number of P-space vectors to add
!> \param offsets specifies the start point in indices and coefficients that defines each vector.
!> \param indices Index in the full space of a contribution to a new P vector
!> \param coefficients Value of a contribution to a new P vector
!> \param pp The P-P block of the matrix, dimensioned (number of existing P + nP, nP)
!> \param parameters On input, the current solution or expansion vector. On exit, the interpolated solution vector.
!> \param action On input, the residual for parameters (non-linear), or action of matrix on parameters (linear). On exit, the expected (non-linear) or actual (linear) residual of the interpolated parameters.
!> \param parametersP On exit, the interpolated solution projected onto the P space.
 SUBROUTINE Iterative_Solver_Add_P(nP,offsets,indices,coefficients,pp,parameters,action,parametersP)
  INTEGER, INTENT(in) :: nP
  INTEGER, INTENT(in), DIMENSION(0:nP) :: offsets
  INTEGER, INTENT(in), DIMENSION(nP) :: indices
  DOUBLE PRECISION, DIMENSION(nP), INTENT(in) :: coefficients
  DOUBLE PRECISION, DIMENSION(*), INTENT(in) :: pp
  DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: parameters
  DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: action
  DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: parametersP
  INTERFACE
   SUBROUTINE IterativeSolverAddPC(nP,offsets,indices,coefficients,pp,parameters,action,parametersP) &
        BIND(C,name='IterativeSolverAddP')
    USE iso_c_binding
    INTEGER(c_size_t), INTENT(in), VALUE :: nP
    INTEGER(c_size_t), INTENT(in), DIMENSION(0:nP) :: offsets
    INTEGER(c_size_t), INTENT(in), DIMENSION(*) :: indices
    REAL(c_double), DIMENSION(*), INTENT(in) :: coefficients
    REAL(c_double), DIMENSION(*), INTENT(in) :: pp
    REAL(c_double), DIMENSION(*), INTENT(inout) :: parameters
    REAL(c_double), DIMENSION(*), INTENT(inout) :: action
    REAL(c_double), DIMENSION(*), INTENT(inout) :: parametersP
   END SUBROUTINE IterativeSolverAddPC
  END INTERFACE
  INTEGER(c_size_t), DIMENSION(0:nP) :: offsetsC
  INTEGER(c_size_t), DIMENSION(SIZE(indices)) :: indicesC
  offsetsC = INT(offsets,c_size_t)
  indicesC = INT(indices,c_size_t)
  CALL IterativeSolverAddPC(INT(nP,c_size_t),offsetsC,indicesC,coefficients, &
   pp,parameters,action,parametersP)
 END SUBROUTINE Iterative_Solver_Add_P

!> \brief the lowest eigenvalues of the reduced problem, for the number of roots sought.
 FUNCTION Iterative_Solver_Eigenvalues()
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Iterative_Solver_Eigenvalues
  INTERFACE
   SUBROUTINE IterativeSolverEigenvalues(eigenvalues) BIND(C,name='IterativeSolverEigenvalues')
    USE iso_c_binding
    REAL(C_double), DIMENSION(*), INTENT(inout) :: eigenvalues
   END SUBROUTINE IterativeSolverEigenvalues
  END INTERFACE
  ALLOCATE (Iterative_Solver_Eigenvalues(m_nroot))
  CALL IterativeSolverEigenvalues(Iterative_Solver_Eigenvalues)
 END FUNCTION Iterative_Solver_Eigenvalues
 
!!> Unit testing of IterativeSolver Fortran binding
 SUBROUTINE Iterative_Solver_Test() BIND(C,name='IterativeSolverFTest')
  INTEGER, PARAMETER :: n=100, nroot=2
  DOUBLE PRECISION, DIMENSION (n,n) :: m
  DOUBLE PRECISION, DIMENSION (n,nroot) :: c,g
  DOUBLE PRECISION, DIMENSION (0,nroot) :: p
  DOUBLE PRECISION, DIMENSION (nroot) :: e,error
  INTEGER :: i,j,root
  PRINT *, 'Test Fortran binding of IterativeSolver'
  m=1
  DO i=1,n
   m(i,i)=3*i
  END DO
  CALL Iterative_Solver_Linear_Eigensystem_Initialize(n,nroot,thresh=1d-8)
  c=0; DO i=1,nroot; c(i,i)=1; ENDDO
  DO i=1,n
   g = MATMUL(m,c)
   CALL Iterative_Solver_Add_Vector(c,g,p)
   e = Iterative_Solver_Eigenvalues()
   WRITE (6,*) 'eigenvalue',e
   DO root=1,nroot
    DO j=1,n
     c(j,root) = c(j,root) - g(j,root)/(m(j,j)-e(i)+1e-15)
    END DO
   END DO
   IF ( Iterative_Solver_End_Iteration(c,g,error)) EXIT
   WRITE (6,*) 'error ',error
  END DO
  CALL Iterative_Solver_Finalize
 END SUBROUTINE Iterative_Solver_Test
END MODULE Iterative_Solver
