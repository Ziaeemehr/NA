program ising2D
!        Metropolis algorithm for the Ising model on a square lattice
IMPLICIT NONE
INTEGER pass,mcs,nequil,N,L,tsave,spin(32,32)
DOUBLE PRECISION w(-8:8),Ce(0:20),Cm(0:20),esave(100),msave(100)
DOUBLE PRECISION E,M,T,accept,cum(10)

L = 1.D0
T = 1.D0
N = L*L

CALL initial(L,T,E,M,w,spin)


! ---------------------------------------------------------
SUBROUTINE initial(L,T,E,M,w,spin)
IMPLICIT NONE


!        seed must not equal 0 to initialize rnd
      dummy = rnd(seed)

!        random initial configuration
DO y = 1,L
    DO x = 1,L       
       IF (rnd(0) .LT. 0.5) THEN
          spin(x,y) = 1
       ELSE
          spin(x,y) = -1
       END IF
          M = M + spin(x,y)
    END DO
END DO

!        compute initial energy
DO y = 1,L
!        periodic boundary conditions
   IF (y .EQ. L) THEN
      up = 1
   ELSE
      up = y + 1
   END IF
   DO x = 1,L
      IF (x .EQ. L) THEN
         right = 1
      ELSE
         right = x + 1
      END IF
      summ = spin(x,up) + spin(right,y)
      E = E - spin(x,y)*summ
   END DO
END DO
!        compute Boltzmann probability ratios
DO dE = -8,8,4
    w(dE) = exp(-dE/T)
END DO
!tsave = 10
END SUBROUTINE

! ---------------------------------------------------------
SUBROUTINE Metropolis(N,L,spin,E,M,w,accept)
!        one Monte Carlo step per spin
IMPLICIT NONE
INTEGER ispin,N,L,x,y,dE,DeltaE,spin(32,32)
DOUBLE PRECISION w(-8:8),E,M,accept
REAL rnd
DO ispin = 1,N
!        random x and y coordinates for trial spin
    x = int(L*rnd(0)) + 1
    y = int(L*rnd(0)) + 1
    dE = DeltaE(x,y,L,spin)
    IF (rnd(0) .LE. w(dE)) THEN
        spin(x,y) = -spin(x,y)
        accept = accept + 1
        M = M + 2*spin(x,y)
        E = E + dE
    END IF
END DO
END SUBROUTINE
! ---------------------------------------------------------
INTEGER FUNCTION DeltaE(x,y,L,spin)
!        periodic boundary conditions
IMPLICIT NONE
INTEGER L,x,y,left,right,up,down,spin(32,32)
IF (x .EQ. 1) THEN
   left = spin(L,y)
   right = spin(2,y)
ELSE IF (x .EQ. L) THEN
   left = spin(L-1,y)
   right = spin(1,y)
ELSE
   left = spin(x-1,y)
   right = spin(x+1,y)
END IF
IF (y .EQ. 1) THEN
   up = spin(x,2)
   down = spin(x,L)
ELSE IF (y .EQ. L) THEN
   up = spin(x,1)
   down = spin(x,L-1)
ELSE
   up = spin(x,y+1)
   down = spin(x,y-1)
END IF
DeltaE = 2*spin(x,y)*(left + right + up + down)
END FUNCTION

SUBROUTINE data(E,M,cum)
!        accumulate data after every Monte Carlo step per spin
IMPLICIT NONE
 DOUBLE PRECISION E,M,cum(10)
 cum(1) = cum(1) + E
 cum(2) = cum(2) + E*E
 cum(3) = cum(3) + M
 cum(4) = cum(4) + M*M
 cum(5) = cum(5) + abs(M)
END SUBROUTINE
! ---------------------------------------------------------
REAL FUNCTION rnd(seed)
!     linear congruential random number generator with shuffling
!     based on ran1 in "Numerical Recipes" second edition
      IMPLICIT NONE
      INTEGER a,m,q,p,n,ndiv,j,k,seed
      REAL rm,rmax
!     m = 2**31 - 1  and  m = a*q + p
      PARAMETER (a = 16807, m = 2147483647, rm = 1.0/m)
      PARAMETER (q = 127773, p = 2836, n = 32, ndiv = 1 + (m-1)/n)
      PARAMETER (rmax = 1.0 - 1.2e-7)
      INTEGER r(n),r0,r1
      SAVE r,r0,r1
      DATA r/n*0/
!     initialize table of random numbers
      IF(seed .NE. 0) THEN
        r1 = abs(seed)
           DO j = n+8,1,-1
             k = r1/q
             r1 = a*(r1-k*q) - p*k
             IF (r1 .LT. 0) r1 = r1 + m
             IF (j .LE. n) r(j) = r1
           ENDDO
        r0 = r(1)
      END IF
!     beginning when not initializing
!     compute r1 = mod(a*r1,m) without overflows
      k = r1/q
      r1 = a*(r1 - k*q) - p*k
      IF (r1 .LT. 0) r1 = r1 + m
      j = 1 + r0/ndiv
      r0 = r(j)
      r(j) = r1
      rnd = min(rm*r0,rmax)
      END FUNCTION
! ---------------------------------------------------------






