PROGRAM LinearEigenSystemExample
 USE IterativeSolverF
 INTEGER, PARAMETER :: n=1000, nroot=2
 DOUBLE PRECISION, DIMENSION (n,n) :: m
 DOUBLE PRECISION, DIMENSION (n) :: c,g
 DOUBLE PRECISION, DIMENSION (nroot) :: e,error
 INTEGER :: i,j
 LOGICAL :: converged
 PRINT *, 'Fortran binding of IterativeSolver'
 m=1; DO i=1,n; m(i,i)=i; END DO
 CALL IterativeSolverLinearEigensystemInitialize(n,nroot)
 c=0; c(1)=1
 DO i=1,n
  g = MATMUL(m,c)
  CALL IterativeSolverLinearEigensystemAddVector(c,g,e)
  DO j=1,n; c(j) = c(j) - g(j)/(m(j,j)-e(1)+1e-15); END DO
  converged = IterativeSolverLinearEigensystemEndIteration(c,g,error)
  PRINT *, 'error =',error,' eigenvalue =',e
  IF (converged) EXIT
 END DO
END PROGRAM LinearEigenSystemExample
