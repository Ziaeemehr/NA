program surface
!surface growth in 1 dimension using dgtsvand Herman_Morison algorithm
Implicit none
INTEGER :: N
INTEGER :: p,iter,nt_step
REAL(8) :: w ! the thickness
REAL(8) , parameter :: nu = 1.d0
REAL(8) , parameter :: D  = 1.d-1
REAL(8),allocatable :: h_old(:),h_new(:), b (:), u(:),v(:), y(:),q(:)
REAL(8),allocatable :: dd(:),du(:),dl(:)
REAL(8) :: dt,dx, alph,L,beta,h_ave
REAL    :: gasdev
INTEGER :: i, j,info,seed

nt_step = 10000
L  = 10.d0
dx = 1.d-2
alph = 5.d-1
! alph = nu * dt / (dx * dx)
dt = alph*dx*dx/nu
print *, alph,dx,dt
N  = L/dx+1
print *, N
allocate(h_old(N),h_new(N),b(N),u(N),v(N),y(N)&
         ,q(N),dd(N),du(N),dl(N))

seed = 12345

open (unit= 21, file = 'wt.txt')
open (unit= 22, file = 'ht.txt')

beta = (2.d0*alph+1.d0)
! define the tridiagonal matrix
dl = -alph
du = -alph
dd = beta
dd(1) = 2.d0 * beta
dd(N) = beta + alph * alph/beta

! Make vector u and v for Sherman Morison algorithm
u = 0.d0
u(1) = -beta
u(N) = -alph
v = 0.d0
v(1) = 1.d0
v(N) = alph/beta

h_old = 0.d0

do iter = 1, nt_step
        do  p = 1, N
            b(p) = sqrt(2 * D * dt) * gasdev(seed) + h_old(p)
        enddo
        call dgtsv(N, 1, dl, dd, du, b, N, info)
        y = b
        ! define the matrix again
        dl = -alph;              du = -alph;         dd = beta
        dd(1) = 2.d0 * beta;     dd(N) = beta + alph * alph/beta
        
        call dgtsv(N, 1, dl, dd, du, u, N, info)
        q = u
        h_new = y - (dot_product(v,y)/(1.d0 + dot_product(v,q))) * q
        h_old = h_new
        ! define the matrix and vector u
        dl = -alph;              du = -alph;         dd = beta
        dd(1) = 2.d0 * beta;     dd(N) = beta + alph * alph/beta	
        u = 0.d0;                u(1) = -beta;       u(N) = -alph
        h_ave = 0.d0
        w    = 0.d0 	
        do i = 1, N
           h_ave = h_ave + h_new(i)
        enddo
        h_ave = h_ave/N
        do i = 1, N
            w = w + (h_new(i)-h_ave)*(h_new(i)-h_ave)
        enddo
        w = sqrt(1/(N*dx)*w)
	write(21,'(es24.15,es24.15)') iter * dt, w
	write(22,'(es24.15,es24.15)') iter * dx, h_ave
enddo

deallocate(h_old,h_new,b,u,v,y,q,dd,du,dl)
end program surface

! -----------------------------------------------
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
  END
! -----------------------------------------------
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
  END





