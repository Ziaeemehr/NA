program surface

Implicit none

INTEGER, parameter :: N = 300
INTEGER :: p,s,it
real(8) :: A(N,N)
real(8) , parameter :: nu = 1.d0
real(8) , parameter :: D  = 1.d-1
real(8) :: h (N), b (N),hh(N),w
real(8) :: dt , dx, Leng, tim, alp, beta(N), X, Y,hbar
real    :: gasdev
integer :: i, j,info, lda, ldb, nrhs,seed
integer, dimension(n) :: ipiv
INTEGER, DIMENSION (8) :: T
nrhs = 1 ! number of right hand sides in b
lda = n  ! leading dimension of a
ldb = n  ! leading dimension of b
dt = 1.d-5
dx = 1.d-1
seed=12345
alp = nu * dt / (dx * dx)

open (unit= 21, file = 'wt.txt')
open (unit= 22, file = 'ht.txt')

h(:) = 0.d0
do s = 1, 1000
    A(:,:) = 0.d0
    do i = 1, N
        do j = 1, N
            if (i==j) A(i,j) = -(2*alp+1)
            if (i==j+1 ) A(i,j) = alp
    	    if (j==i+1 ) A(i,j) = alp
        enddo
    enddo
    A(N,1) = alp
    A(1,N) = alp
        do  p = 1, N
	    b(p) = -sqrt(2 * D * dt) * gasdev(seed) - h(p)
	enddo
	call dgesv(N, nrhs, A, lda, ipiv, b, ldb, info)
	h(1:N) = b(1:N)
        hbar = 0.d0
	w    = 0.d0
        !hbar = sum(h(1:N))
	!hbar = hbar/N    
	!hh(1:N) = h(1:N)-hbar
	!w = sqrt(1/(N*dx)*sum(hh(1:N)*hh(1:N)))
        do i = 1, N
	   hbar = hbar + h(i)
        enddo
        hbar = hbar/N
        do i = 1, N
            w = w + (h(i)-hbar)*(h(i)-hbar)
        enddo
        w = sqrt(1.d0/(N*dx)*w)
	write(21,'(es24.15,es24.15)') s*dt, w
	write(22,'(es24.15,es24.15)') s*dx, hbar
enddo


end program surface

!==========================================================
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





