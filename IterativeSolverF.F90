!> @brief IterativeSolver Fortran binding
MODULE IterativeSolverF
!> @example LinearEigensystemExampleF.F90
 USE iso_c_binding
 !PRIVATE
 INTEGER(c_size_t) :: m_nq, m_np, m_nroot

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

 SUBROUTINE IterativeSolverLinearEigensystemInitialize(nq,nroot)
  INTEGER, INTENT(in) :: nq
  INTEGER, INTENT(in) :: nroot
  INTERFACE
   SUBROUTINE IterativeSolverLinearEigensystemInitializeC(nq,nroot) &
     BIND(C,name='IterativeSolverLinearEigensystemInitialize')
    USE iso_c_binding
    INTEGER(C_size_t), INTENT(in), VALUE :: nq
    INTEGER(C_size_t), INTENT(in), VALUE :: nroot
   END SUBROUTINE IterativeSolverLinearEigensystemInitializeC
  END INTERFACE
  m_nq=int(nq,kind=c_size_t)
  m_nroot=INT(nroot,kind=c_size_t)
  CALL IterativeSolverLinearEigensystemInitializeC(m_nq,m_nroot)
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
 INTEGER, PARAMETER :: n=100, nroot=1
 DOUBLE PRECISION, DIMENSION (n,n) :: m
 DOUBLE PRECISION, DIMENSION (n) :: c,g
 DOUBLE PRECISION, DIMENSION (nroot) :: e,error
 INTEGER :: i,j
 PRINT *, 'Test Fortran binding of IterativeSolver'
 m=1
 DO i=1,n
  m(i,i)=i
 END DO
 CALL IterativeSolverLinearEigensystemInitialize(n,1)
 c=0
 c(1)=1
 DO i=1,n
  g = MATMUL(m,c)
  CALL IterativeSolverLinearEigensystemAddVector(c,g,e)
  WRITE (6,*) 'eigenvalue',e
  DO j=1,n
   c(j) = c(j) - g(j)/(m(j,j)-e(1)+1e-15)
  END DO
  IF ( IterativeSolverLinearEigensystemEndIteration(c,g,error)) EXIT
  WRITE (6,*) 'error ',error
 END DO
 END SUBROUTINE IterativeSolverFTest
END MODULE IterativeSolverF
