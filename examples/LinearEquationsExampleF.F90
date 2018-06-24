!> @example LinearEquationsExampleF.F90
!> This is an example of simplest use of the LinearEquations framework for iterative
!> finding of the lowest few eigensolutions of a large matrix.
PROGRAM LinearEigenSystemExample
 USE Iterative_Solver
 INTEGER, PARAMETER :: n=300, nroot=2
 DOUBLE PRECISION, DIMENSION (n,n) :: m
 DOUBLE PRECISION, DIMENSION (n,nroot) :: c,g, rhs
 DOUBLE PRECISION, DIMENSION (nroot) :: e,error
 DOUBLE PRECISION, PARAMETER :: alpha=300
 DOUBLE PRECISION :: aughes
 DOUBLE PRECISION, DIMENSION(5), PARAMETER :: augmented_hessian_factors = [0.0, .001, .01, .1, 1.0]
 INTEGER :: i,j,root
 LOGICAL :: converged
 PRINT *, 'Fortran binding of IterativeSolver'
 DO i=1,n; m(i,i)=alpha*i+2*i-2; DO j=1,n; if (i.ne.j) m(i,j)=i+j-2; END DO; END DO
 do i=1,nroot
 do j=1,n
 rhs(j,i) = 1/dble(j+i-1)
 end do
 end do
 do iaug=1,size(augmented_hessian_factors)
 aughes = augmented_hessian_factors(iaug)
 print *, 'solve linear system with augmented hessian factor ',aughes
 CALL Iterative_Solver_Linear_Equations_Initialize(n,nroot,rhs,aughes,thresh=1d-11,verbosity=1)
 c=0; DO i=1,nroot; c(i,i)=1; ENDDO
 DO i=1,n
  g = MATMUL(m,c)
  CALL Iterative_Solver_Add_Vector(c,g,e)
  DO root=1,nroot
   DO j=1,n
    c(j,root) = c(j,root) - g(j,root)/(m(j,j))
   END DO
  END DO
  converged = Iterative_Solver_End_Iteration(c,g,error)
  IF (converged) EXIT
 END DO
 PRINT *, 'error =',error
 do i=1,nroot
 print *, 'solution ',c(1:min(n,10),i)
 end do
 CALL Iterative_Solver_Finalize
 enddo
END PROGRAM LinearEigenSystemExample
