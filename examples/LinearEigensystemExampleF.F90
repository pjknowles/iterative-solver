PROGRAM LinearEigenSystemExample
 USE IterativeSolverF
 INTEGER, PARAMETER :: n=100, nroot=1
 DOUBLE PRECISION, DIMENSION (n,n) :: m
 DOUBLE PRECISION, DIMENSION (n) :: c,g
 DOUBLE PRECISION, DIMENSION (nroot) :: e,error
 INTEGER :: i,j
 PRINT *, 'Fortran binding of IterativeSolver'
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
END PROGRAM LinearEigenSystemExample
