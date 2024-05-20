! Created by Peter Knowles on 20/05/2024.

program LinearEquationsSmall
  use Iterative_Solver
  use iso_fortran_env, only: dp => real64
  implicit double precision (a-h,o-z)
  integer, parameter :: iout=6
  real(dp), pointer, contiguous :: a(:,:),rhs(:),b(:),at(:,:),ab(:)
  real(dp), pointer, contiguous :: error(:),tmp1(:,:),tmp2(:,:)
  integer, pointer, contiguous :: ipiv(:)

  n=10
  nrhs=1

  allocate(a(n,n))
  allocate(at(n,n))
  allocate(ab(n))
  allocate(b(n))
  allocate(rhs(n))
  allocate(ipiv(n))
  allocate(error(nrhs))
  allocate(tmp1(n,nrhs))
  allocate(tmp2(n,nrhs))

  call random_seed()
  call random_number(a)
  call random_number(rhs)
  do i=1,n
    do j=1,n
      a(i,j) = i+j*100
    end do
  end do
  rhs=1

  !make a diagonal dominant
  do i=1,n
    do j=1,i
      if(i.eq.j) then
        a(i,i)=a(i,i)*1d1
      else
        a(i,j)=a(i,j)*1d-1
        a(j,i)=a(j,i)*1d-1
      endif
    enddo
  enddo
  call outsqr(a,n,n,n,'a')
  call outvec(rhs,n,'rhs')

  !use standard solver for x.A=b
  b(:)=rhs(:)
  at(:,:)=transpose(a)
  call dgesv(n,nrhs,at,n,ipiv,b,n,info)
  call outvec(b,n,'solution b from dgesv')

  !iterative solver
  call Iterative_Solver_Linear_Equations_Initialize(n,nrhs,rhs,0d0, &
      & hermitian=.false.,thresh=1d-12,verbosity=1)

  !starting guess
  b(:)=0d0
  b(1)=1d0

  write(iout,'(/,t5,a,t20,a)') 'iteration','error'
  do iter=1,n
    ab=matmul(b,a) !b.A
    !call outvec(b,n,'b')
    !call outvec(ab,n,'ab')

    tmp1(:,1)=b
    tmp2(:,1)=ab
    nwork=Iterative_Solver_Add_Vector(tmp1,tmp2,.false.)
    b=tmp1(:,1)
    ab=tmp2(:,1)

    do i=1,n
      if(abs(a(i,i)).gt.1d-6) then
        b(i)=b(i)-ab(i)/a(i,i)
      endif
    enddo

    tmp1(:,1)=b
    tmp2(:,1)=ab
    nwork=Iterative_Solver_End_Iteration(tmp1,tmp2)
    error=Iterative_Solver_Errors()
    b=tmp1(:,1)
    ab=tmp2(:,1)
    errmax=maxval(abs(error))
    write(iout,'(t5,i5,t15,e14.5,t32,i4)') iter,errmax
    if(nwork.le.0) exit
  enddo

  call outvec(b,n,'solution b from IterativeSolver')

  call Iterative_Solver_Finalize()




end program LinearEquationsSmall

subroutine outsqr (q,idim,ia,ib,title)
  implicit double precision(a-h,o-z)
  integer, parameter :: iout=6
  character(len=*) title
  dimension q(idim,*)
  double precision x(8)
  write (iout,'(/1X,(A))') title
  nc=8
  m=1
  n=nc
  10    if (ib.lt.m) return
  n=min0(n,ib)
  write (iout,30) (i,i=m,n)
  do 20 j=1,ia
    ii=0
    do i=m,n
      ii=ii+1
      x(ii)=q(j,i)
      if(abs(x(ii)).lt.1.d-8) x(ii)=0  !remove negative zeros from output disturb in diffs)
    end do
    write (iout,50) j,(x(i) ,i=1,ii)
  20    continue
  m=m+nc
  n=n+nc
  goto 10
  30    format(1x,99i14)
  40    format(/)
  50    format(2x,i3,99f14.8)
end
subroutine outvec (p,ib,title)
  implicit double precision(a-h,o-z)
  integer, parameter :: iout=6
  character(len=*) title
  dimension p(*),x(10)
  write (iout,'(/1X,(A))') title
  m=1
  n=10
  10    if (ib.lt.m) return
  n=min0(n,ib)
  ii=0
  do i=m,n
    ii=ii+1
    x(ii)=p(i)
    if(abs(x(ii)).lt.1.d-8) x(ii)=0  !remove negative zeros from output disturb in diffs)
  end do
  if(ib.lt.1000) then
    write(iout,21) m,n,(x(i) ,i=1,ii)
  else
    write(iout,20) m,n,(x(i) ,i=1,ii)
  endif
  m=m+10
  n=n+10
  goto 10
  20    format(1x,i0,'-',i0,10f12.6)
  21    format(1x,i5,'-',i4,10f12.6)
end