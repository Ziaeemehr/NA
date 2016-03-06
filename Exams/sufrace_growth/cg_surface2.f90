program surface
!surface growth in 1 dimension using Conjugate Gradient Algorithm
Implicit none
INTEGER,PARAMETER :: N = 100
INTEGER :: p,iter,nt_step
REAL(8) :: w ! thickness
REAL(8) , parameter :: nu = 1.d0
REAL(8) , parameter :: D  = 1.d-1
REAL(8) :: h_old(N),h_new(N), b (N),h(N),g_new(N),g(N)
REAL(8) :: dt,dx, alph,L,beta,h_ave,tol,s,a
REAL(8) :: denom,gg_new,gg,gnorm
REAL    :: gasdev
INTEGER :: i, j,info,seed

L  = 2.d0
dx = L/dble(N-1)
a = 5.d-1
dt = a*dx*dx/nu
tol = 1.d-3
print *, a,dx,dt
! stop


seed = 12345
open (unit= 21, file = 'wt.txt')
open (unit= 22, file = 'ht.txt')

s = (2.d0 * a + 1.d0)
h_old = 0.d0
! DO i = 1,N
!    h_new(i) = - sqrt(2*D*dt)*gasdev(seed)
! ENDDO
! print *, h_new
h_new = 0.

DO iter = 2, 10
    g(1) = s*h_new(1)-a*h_new (2)-a*h_new(N)-h_old(1) &
          - sqrt(2*D*dt)*gasdev(seed)
    DO i = 2,N-1
        g(i) = - a*h_new(i-1)+s*h_new(i)-a*h_new(i+1)-h_old(i)&
          - sqrt(2*D*dt)*gasdev(seed)
    ENDDO
    g(N) = -a*h_new(1)-a*h_new(N-1)+s*h_new(N)-h_old(N)&
          - sqrt(2*D*dt)*gasdev(seed)
! -----------------------------------------------    
    ! CG LOOP    
    h = -g
!     print*, h
!     stop
    DO  p = 1,200
       denom = -2.d0*a*h(1)*h(N)+s*h(N)*h(N)
       DO i = 1,N-1
           denom = denom+s*h(i)*h(i)-2.d0*a*h(i)*h(i+1)
       ENDDO
       gg = DOT_PRODUCT(g,g)
       alph = gg / denom
       h_new=h_new+alph*h
! -----------------------------------------------
       ! update g
       g_new(1) = g(1) + alph*(s*h(1) - a*h(2) - a*h(N))
       DO i = 2,N-1
           g_new(i) = g(i) + alph * (-a*h(i-1)-a*h(i+1)+s*h(i))
       ENDDO
       g_new(N) = g(N) + alph * (-a *h(1) - a* h(N-1) * s * h(N))
! -----------------------------------------------       
       gg_new = DOT_PRODUCT(g_new,g_new)
       beta = gg_new/gg
       h = -g_new + beta * h
       gnorm = SQRT(gg_new)
       print *, p,gnorm
       IF (gnorm < tol) EXIT
       g = g_new
    ENDDO
! -----------------------------------------------    
    h_old = h_new
    h_ave = 0.d0
    w = 0.d0
    DO i = 1,N
        h_ave = h_ave + h_new(i)
    ENDDO
    h_ave = h_ave / N    
    DO i = 1,N
        w = w + (h_new(i)-h_ave)*(h_new(i)-h_ave)
    ENDDO
    w = SQRT(1.d0/(N*dx)*w)
    write(21,'(es24.15,es24.15)') iter * dt, w
    write(22,'(es24.15,es24.15)') iter * dx, h_ave
ENDDO

END PROGRAM surface

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
