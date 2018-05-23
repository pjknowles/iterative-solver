PROGRAM LinearEigenSystemExample
 USE IterativeSolverF
 INTEGER, PARAMETER :: n=1000, nroot=2
 DOUBLE PRECISION, DIMENSION (n,n) :: m
 DOUBLE PRECISION, DIMENSION (n,nroot) :: c,g
 DOUBLE PRECISION, DIMENSION (nroot) :: e,error
 INTEGER :: i,j,root
 LOGICAL :: converged
 PRINT *, 'Fortran binding of IterativeSolver'
 m=1; DO i=1,n; m(i,i)=i; END DO
 CALL IterativeSolverLinearEigensystemInitialize(n,nroot)
 c=0; DO i=1,nroot; c(i,i)=1; ENDDO
 DO i=1,n
  g = MATMUL(m,c)
  CALL IterativeSolverLinearEigensystemAddVector(c,g,e)
  DO root=1,nroot
   DO j=1,n; c(j,root) = c(j,root) - g(j,root)/(m(j,j)-e(root)+1e-15); END DO
  END DO
  converged = IterativeSolverLinearEigensystemEndIteration(c,g,error)
  PRINT *, 'error =',error,' eigenvalue =',e
  IF (converged) EXIT
 END DO
END PROGRAM LinearEigenSystemExample
