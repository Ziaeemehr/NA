program landao
!
! To solve the Landao- Ginzberg eauation for a 1 dimensional loop
! m(x,t) = r and periodic boundary condition 
! Variables :
! dl, dd, du subdiagonal, maindiagonal and superdiagonal elements respectively
! To plot the result:
! echo 'splot "xt.txt" nonuniform matrix w l lt 7' | gnuplot --persist
! echo 'plot "x.txt" u 1:2 w l lt 7' | gnuplot --persist
! -------------------------------------------------------------------

implicit none
integer :: n ,m
integer :: i,j,seed,info
real(8) :: alph, l, beta
real(8) :: d_x,d_t,K,a
real(8),allocatable :: m_new(:),m_old(:),dl(:),dd(:),du(:),y(:)
real(8),allocatable :: u(:),v(:),q(:)
real    :: gasdev

l = 3.      ! change l over [2,3] to present the phase transition
!time = 1.d-1
!d_x = l/real(n-1)
d_x = 5.d-2
n = l/d_x + 1
m = 100
!d_t  = time/real(m-1)
d_t = 125.d-5
write(*, 99)  "n= ",n," m= ",m,' dx= ',d_x,'dt= ',d_t
99 format (A5,I3,A5,I3,A5,f7.4,A5,f10.7)

allocate (m_new(n),m_old(n),dl(n),du(n),dd(n),y(n),u(n),v(n),q(n))



K = 1.d0
a = 1.d0
alph = (K * d_t)/(d_x*d_x)
print *, 'alpha', alph

if (alph > 0.5) then
print*, "Solution may be unaccurate"
!stop
endif

! Define the tridiagonal matrix A'
beta = 1.d0 + 2.d0 * alph + a * d_t
du = -alph
dl = -alph
dd = beta
dd(1) = 2.d0 * beta
dd(n) = beta + alph*alph/beta

! But the matrix is tridiagonal + 2 extra elements
! A(1,N) = A(N,1) = -alph
seed = 12345
do i = 1,n 
m_old(i) = gasdev(seed)
enddo

! Define tow Morison vectors
u = 0.d0
u(1) = - beta
u(n) = -alph
v = 0.d0
v(1) = 1.d0
v(n) = alph/beta


! print the x positions at the first row
open(unit=7, file="xt.txt")
open(unit=14, file="x.txt")
write(7,120) n,((i-1)*d_x,i=1,n)
write(7,110) 0.d0,(m_old(i),i=1,n)
110 format (f6.3,999f12.6)
120 format (I3,999f12.2)

! main loop for each time step
do j = 2,m-1   
   call dgtsv( n , 1 , dl , dd , du , m_old , n , info )
   y = m_old
   ! define the matrix A' again
   du = -alph;   dl = -alph;   dd = beta;   dd(1) = 2.d0 * beta;   dd(n) = beta + alph*alph/beta
   call dgtsv( n , 1 , dl , dd , du , u , n , info )
   q = u
   m_new = y - dot_product(v,y)/(1+dot_product(v,q)) * q
   ! define the matrix A' again and u
   du = -alph; dl = -alph;    dd = beta;   dd(1) = 2.d0 * beta;   dd(n) = beta + alph*alph/beta   
   u = 0.d0;   u(1) = - beta; u(n) = -alph
   if (j==60) then
       do i=1,n
           write(14,'(f10.4,f14.7)'),  (i-1)*d_x,m_new(i)
       enddo
   endif
   write(7,110) j*d_t ,(m_new(i),i=1,n)
   m_old = m_new
enddo

! Printing the results
!write(*,120) ((i-1)*h,i=2,nx-1)
!120 format (7x,' t',7x,'x = ',f5.2,9f11.2)


close(7)

end program
! -------------------------------------------------------------------

  FUNCTION gasdev(idum)
  INTEGER idum
  REAL gasdev
    !CU    USES ran1
  INTEGER iset
  REAL fac,gset,rsq,v1,v2,ran1
  SAVE iset,gset
  DATA iset/0/
  if (idum.lt.0) iset=0
  if (iset.eq.0) then
1       v1=2.*ran1(idum)-1.
    v2=2.*ran1(idum)-1.
    rsq=v1**2+v2**2
    if(rsq.ge.1..or.rsq.eq.0.)goto 1
    fac=sqrt(-2.*log(rsq)/rsq)
    gset=v1*fac
    gasdev=v2*fac
    iset=1
  else
    gasdev=gset
    iset=0
  endif
  return
  END FUNCTION
! -------------------------------------------------------------------

  FUNCTION ran1(idum)
  INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
  REAL ran1,AM,EPS,RNMX
  PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
  NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
  INTEGER j,k,iv(NTAB),iy
  SAVE iv,iy
  DATA iv /NTAB*0/, iy /0/
  if (idum.le.0.or.iy.eq.0) then
    idum=max(-idum,1)
    do 11 j=NTAB+8,1,-1
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      if (j.le.NTAB) iv(j)=idum
11      continue
    iy=iv(1)
  endif
  k=idum/IQ
  idum=IA*(idum-k*IQ)-IR*k
  if (idum.lt.0) idum=idum+IM
  j=1+iy/NDIV
  iy=iv(j)
  iv(j)=idum
  ran1=min(AM*iy,RNMX)
  return
  END FUNCTION
  

