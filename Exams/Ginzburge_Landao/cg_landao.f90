PROGRAM LANDAO_CONJUGATE_GRADIENT
!
! To solve the Landao- Ginzberg eauation for a 1 dimensional loop
! m(x,t) = r and periodic boundary condition 
! To plot the result:
! echo 'splot "xt.txt" nonuniform matrix w l lt 7' | gnuplot --persist
! echo 'plot "x.txt" u 1:2 w l lt 7' | gnuplot --persist
! -------------------------------------------------------------------

implicit none
integer :: N ,M
integer :: i,j,seed,info,iter,ip
real(8) :: alph,l,beta,s,a,K,aa,gg,denom,gnorm,tol
real(8) :: dx,dt
real(8),allocatable :: m_new(:),m_old(:),g(:),g_new(:),h(:)
real    :: gasdev

tol = 1.d-6
l = 2.7      ! change l over [2,3] to present the phase transition
dx = 5.d-2
N = l/dx + 1
m = 10
K = 1.d0
a = 1.d0
aa = 5.d-1
dt = aa * dx*dx/K
write(*, 99)  "n= ",n," m= ",m,' dx= ',dx,'dt= ',dt
99 format (A5,I3,A5,I3,A5,f7.4,A5,f10.7)

allocate (m_new(n),m_old(n),h(N),g_new(N),g(N))


s = 1.d0 + 2.d0 * aa + a * dt
seed = 12345
do i = 1,N
    m_old(i) = gasdev(seed)    
enddo

! print the x positions at the first row
open(unit=7, file="xt.txt")
open(unit=14, file="x.txt")

write(7,110) N,((i-1)*dx,i=1,N)
write(7,120) 0.d0,(m_old(i),i=1,N)
110 format (I12,*(f12.6))
120 format (f12.6,*(f12.6))

m_new = m_old

DO iter = 1,100 ! time iteration
    g(1) = s * m_new(1)- aa*m_new(2)- aa*m_new(N)- m_old(1)
    DO i = 2,N-1
       g(i) = -aa*m_new(i-1)+s*m_new(i)-aa*m_new(i+1)-m_old(i)
    ENDDO   
    g(N) = -aa*m_new(1)-aa*m_new(N-1)+s*m_new(N)-m_old(N)
    ! CG loops
    h = -g
    DO ip = 1,100
        denom = -2*aa*h(1)*h(N) + s*h(N)*h(N)
        do i = 1,N-1
            denom = denom + s*h(i)*h(i)-2.d0*aa*h(i)*h(i+1)
        enddo
        gg = DOT_PRODUCT(g,g)
        alph = gg / denom
        m_new = m_new + alph * h        
        ! update g
        g_new(1) = g(1) + alph * (s*h(1)-aa*h(2)-aa*h(N))
        do i = 2,N-1
           g_new(i) = g(i) + alph * (-aa*h(i-1)+ s*h(i) - aa*h(i+1))
        enddo
        g_new(N) = g(N) + alph * (-aa*h(1)-aa*h(N-1)+s*h(N))
        
        beta = DOT_PRODUCT(g_new,g_new)
        h = -g_new + beta * h
        gnorm = SQRT(DOT_PRODUCT(g_new,g_new))
!         print *, ip , gnorm
        IF (gnorm < tol) EXIT          
        g = g_new
    ENDDO
     write(7,120) iter*dt ,(m_new(i),i=1,N)
     m_old = m_new
     IF (iter == 50) THEN
     DO i= 1,N
         write(14,'(f12.6,f12.6)'),  (i-1)*dx,m_new(i) 
     ENDDO
     ENDIF
ENDDO

close(7 )
close(14)

END PROGRAM
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