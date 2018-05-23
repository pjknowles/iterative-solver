!> @brief IterativeSolver Fortran binding
MODULE IterativeSolverF
 USE iso_c_binding
 PUBLIC :: IterativeSolverLinearEigensystemInitialize
 PUBLIC :: IterativeSolverLinearEigensystemAddVector, IterativeSolverLinearEigensystemEndIteration
 PRIVATE
 INTEGER(c_size_t) :: m_nq, m_nroot

 INTERFACE
  SUBROUTINE IterativeSolverLinearEigensystemAddVector(c,g,e) BIND(C,name='IterativeSolverLinearEigensystemAddVector')
   USE iso_c_binding
   DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: c
   DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: g
   DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: e
  END SUBROUTINE IterativeSolverLinearEigensystemAddVector
 END INTERFACE

 INTERFACE
  LOGICAL FUNCTION IterativeSolverLinearEigensystemEndIteration(c,g,error) &
       BIND(C,name='IterativeSolverLinearEigensystemEndIteration')
   USE iso_c_binding
   DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: c
   DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: g
   DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: error
  END FUNCTION IterativeSolverLinearEigensystemEndIteration
 END INTERFACE

CONTAINS

!> @example LinearEigensystemExampleF.F90
!> \brief Finds the lowest eigensolutions of a matrix using Davidson's method, i.e. preconditioned Lanczos
!>
!> Example of simplest use: @include LinearEigensystemExampleF.F90
 SUBROUTINE IterativeSolverLinearEigensystemInitialize(nq,nroot,thresh,maxIterations,verbosity)
  INTEGER, INTENT(in) :: nq !< dimension of matrix
  INTEGER, INTENT(in) :: nroot !< number of eigensolutions desired
  DOUBLE PRECISION, INTENT(in), OPTIONAL :: thresh !< convergence threshold
  INTEGER, INTENT(in), OPTIONAL :: maxIterations !< maximum number of iterations
  INTEGER, INTENT(in), OPTIONAL :: verbosity !< how much to print. Default is zero, which prints nothing except errors
  INTERFACE
   SUBROUTINE IterativeSolverLinearEigensystemInitializeC(nq,nroot,thresh,maxIterations,verbosity) &
        BIND(C,name='IterativeSolverLinearEigensystemInitialize')
    USE iso_c_binding
    INTEGER(C_size_t), INTENT(in), VALUE :: nq
    INTEGER(C_size_t), INTENT(in), VALUE :: nroot
    DOUBLE PRECISION, INTENT(in), VALUE :: thresh
    INTEGER(C_int), INTENT(in), VALUE :: maxIterations
    INTEGER(C_int), INTENT(in), VALUE :: verbosity
   END SUBROUTINE IterativeSolverLinearEigensystemInitializeC
  END INTERFACE
  INTEGER(c_int) :: verbosityC, maxIterationsC
  DOUBLE PRECISION :: threshC
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
  CALL IterativeSolverLinearEigensystemInitializeC(m_nq,m_nroot,threshC, maxIterationsC, verbosityC)
 END SUBROUTINE IterativeSolverLinearEigensystemInitialize
 
!!> Add P-space vectors to the expansion set
 SUBROUTINE IterativeSolverAddP(indices,coefficients,pp)
  INTEGER, INTENT(in), DIMENSION(:) :: indices
  DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: coefficients
  DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: pp
  INTERFACE
   SUBROUTINE IterativeSolverAddPC(indices,coefficients,pp) &
        BIND(C,name='IterativeSolverAddP')
    USE iso_c_binding
    INTEGER(c_size_t), INTENT(in), DIMENSION(*) :: indices
    DOUBLE PRECISION, DIMENSION(*), INTENT(in) :: coefficients
    DOUBLE PRECISION, DIMENSION(*), INTENT(in) :: pp
   END SUBROUTINE IterativeSolverAddPC
  END INTERFACE
  INTEGER(c_size_t), DIMENSION(SIZE(indices)) :: indicesC
  indicesC = INT(indices,c_size_t)
  CALL IterativeSolverAddPC(indicesC,coefficients,pp)
 END SUBROUTINE IterativeSolverAddP
 
!!> Unit testing of IterativeSolver Fortran binding
 SUBROUTINE IterativeSolverFTest() BIND(C,name='IterativeSolverFTest')
  INTEGER, PARAMETER :: n=100, nroot=2
  DOUBLE PRECISION, DIMENSION (n,n) :: m
  DOUBLE PRECISION, DIMENSION (n,nroot) :: c,g
  DOUBLE PRECISION, DIMENSION (nroot) :: e,error
  INTEGER :: i,j,root
  PRINT *, 'Test Fortran binding of IterativeSolver'
  m=1
  DO i=1,n
   m(i,i)=i
  END DO
  CALL IterativeSolverLinearEigensystemInitialize(n,nroot)
  c=0; DO i=1,nroot; c(i,i)=1; ENDDO
  DO i=1,n
   g = MATMUL(m,c)
   CALL IterativeSolverLinearEigensystemAddVector(c,g,e)
   WRITE (6,*) 'eigenvalue',e
   DO root=1,nroot
    DO j=1,n
     c(j,root) = c(j,root) - g(j,root)/(m(j,j)-e(i)+1e-15)
    END DO
   END DO
   IF ( IterativeSolverLinearEigensystemEndIteration(c,g,error)) EXIT
   WRITE (6,*) 'error ',error
  END DO
 END SUBROUTINE IterativeSolverFTest
END MODULE IterativeSolverF
