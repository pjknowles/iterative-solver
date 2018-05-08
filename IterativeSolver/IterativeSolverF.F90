!> @brief IterativeSolver Fortran binding
MODULE IterativeSolverF
 USE iso_c_binding
 !PRIVATE
 INTEGER(c_size_t) :: m_nq, m_np, m_nroot

 INTERFACE
  SUBROUTINE IterativeSolverDavidsonInterpolate(c,g,e) BIND(C,name='IterativeSolverDavidsonInterpolate')
   USE iso_c_binding
   DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: c
   DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: g
   DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: e
  END SUBROUTINE IterativeSolverDavidsonInterpolate
 END INTERFACE

 INTERFACE
  LOGICAL FUNCTION IterativeSolverDavidsonFinalize(c,g,error) BIND(C,name='IterativeSolverDavidsonFinalize')
   USE iso_c_binding
   DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: c
   DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: g
   DOUBLE PRECISION, DIMENSION(*), INTENT(inout) :: error
  END FUNCTION IterativeSolverDavidsonFinalize
 END INTERFACE

CONTAINS

 SUBROUTINE IterativeSolverDavidsonInitialize(nq,np,nroot,pp)
  INTEGER, INTENT(in) :: nq
  INTEGER, INTENT(in) :: np
  INTEGER, INTENT(in) :: nroot
  DOUBLE PRECISION, DIMENSION(np,*), INTENT(in) :: pp
  INTERFACE
   SUBROUTINE IterativeSolverDavidsonInitializeC(nq,np_,nroot,pp_) BIND(C,name='IterativeSolverDavidsonInitialize')
    USE iso_c_binding
    INTEGER(C_size_t), INTENT(in), VALUE :: nq
    INTEGER(C_size_t), INTENT(in), VALUE :: np_
    INTEGER(C_size_t), INTENT(in), VALUE :: nroot
    DOUBLE PRECISION, DIMENSION(*), INTENT(in) :: pp_
   END SUBROUTINE IterativeSolverDavidsonInitializeC
  END INTERFACE
  m_nq=int(nq,kind=c_size_t)
  m_np=INT(np,kind=c_size_t)
  m_nroot=INT(nroot,kind=c_size_t)
  CALL IterativeSolverDavidsonInitializeC(m_nq,m_np,m_nroot,pp)
 END SUBROUTINE IterativeSolverDavidsonInitialize

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
 CALL IterativeSolverDavidsonInitialize(n,0,1,m)
 c=0
 c(1)=1
 DO i=1,n
  g = MATMUL(m,c)
  CALL IterativeSolverDavidsonInterpolate(c,g,e)
  WRITE (6,*) 'eigenvalue',e
  DO j=1,n
   c(j) = c(j) - g(j)/(m(j,j)-e(1)+1e-15)
  END DO
  IF ( IterativeSolverDavidsonFinalize(c,g,error)) EXIT
  WRITE (6,*) 'error ',error
 END DO
 END SUBROUTINE IterativeSolverFTest
END MODULE IterativeSolverF
