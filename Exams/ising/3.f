      PROGRAM ising
*     Metropolis algorithm for the Ising model on a square lattice
      IMPLICIT NONE
      INTEGER pass,mcs,nequil,N,L,tsave,spin(32,32)
      DOUBLE PRECISION w(-8:8),Ce(0:20),Cm(0:20),esave(100),msave(100)
      DOUBLE PRECISION E,M,T,accept,cum(10)
      
      CALL initial(N,L,T,mcs,nequil,tsave,E,M,w,spin)
*     equilibrate system
      DO pass = 1,nequil
         CALL Metropolis(N,L,spin,E,M,w,accept)
      END DO
      CALL initialize(cum,Ce,Cm,esave,msave,accept)
*     accumulate data while updating spins
      DO pass = 1,mcs
         CALL Metropolis(N,L,spin,E,M,w,accept)
         CALL data(E,M,cum)
         CALL correl(Ce,Cm,E,M,esave,msave,pass,tsave)
      END DO
      CALL output(T,N,mcs,cum,accept)
      CALL c_output(Ce,Cm,cum,mcs,tsave)
      END

      SUBROUTINE initial(N,L,T,mcs,nequil,tsave,E,M,w,spin)
      IMPLICIT NONE
      INTEGER mcs,nequil,tsave,N,L,x,y,up,right,dE,seed,spin(32,32)
      DOUBLE PRECISION w(-8:8),E,M,T,sum
      REAL dummy,rnd
      
      WRITE(*,*)'linear dimension of lattice = '
      READ(*,*) L
      WRITE(*,*)'reduced temperature = '
      READ(*,*) T
      N = L*L
      WRITE(*,*) '# MC steps per spin for equilibrium = '
      READ(*,*) nequil
      WRITE(*,*) '# MC steps per spin = '
      READ(*,*) mcs
      WRITE(*,*) 'random number seed = '
      READ(*,*) seed
*     seed must not equal 0 to initialize rnd
      dummy = rnd(seed)
*     random initial configuration
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
*     compute initial energy
      DO y = 1,L
*        periodic boundary conditions
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
            sum = spin(x,up) + spin(right,y)
            E = E - spin(x,y)*sum
         END DO
      END DO
*     compute Boltzmann probability ratios
      DO dE = -8,8,4
         w(dE) = exp(-dE/T)
      END DO
      tsave = 10
      END

      SUBROUTINE initialize(cum,Ce,Cm,esave,msave,accept)
      IMPLICIT NONE
      INTEGER i
      DOUBLE PRECISION Ce(0:20),Cm(0:20),esave(100),msave(100)
      DOUBLE PRECISION accept,cum(10)
      DO i = 1,100
         esave(i) = 0
         msave(i) = 0
      END DO
      DO i = 1,20
         Ce(i) = 0
         Cm(i) = 0
      END DO
      DO i = 1,5
         cum(i) = 0
      END DO
      accept = 0
      END

      SUBROUTINE Metropolis(N,L,spin,E,M,w,accept)
*     one Monte Carlo step per spin
      IMPLICIT NONE
      INTEGER ispin,N,L,x,y,dE,DeltaE,spin(32,32)
      DOUBLE PRECISION w(-8:8),E,M,accept
      REAL rnd
      DO ispin = 1,N
*        random x and y coordinates for trial spin
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
      END

      INTEGER FUNCTION DeltaE(x,y,L,spin)
*     periodic boundary conditions
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
      END

      SUBROUTINE data(E,M,cum)
*     accumulate data after every Monte Carlo step per spin
      IMPLICIT NONE
      DOUBLE PRECISION E,M,cum(10)
      cum(1) = cum(1) + E
      cum(2) = cum(2) + E*E
      cum(3) = cum(3) + M
      cum(4) = cum(4) + M*M
      cum(5) = cum(5) + abs(M)
      END

      SUBROUTINE output(T,N,mcs,cum,accept)
*     average per spin
      IMPLICIT NONE
      INTEGER N,mcs
      DOUBLE PRECISION T,accept,norm,eave,e2ave,mave,m2ave
      DOUBLE PRECISION abs_mave,cum(10)
      norm = 1.0/(mcs*N)
      accept = accept*norm
      eave = cum(1)*norm
      e2ave = cum(2)*norm
      mave = cum(3)*norm
      m2ave = cum(4)*norm
      abs_mave = cum(5)*norm
      WRITE(*,*)
      WRITE(*,'(A14,F8.4)') 'Temperature = ',T
      WRITE(*,'(A25,F8.4)') 'acceptance probability = ',accept
      WRITE(*,'(A23,F8.4)') 'mean energy per spin = ', eave
      WRITE(*,'(A30,F10.4)') 'mean squared energy per spin = ', e2ave
      WRITE(*,'(A29,F8.4)') 'mean magnetization per spin = ', mave
      WRITE(*,'(A25,F8.4)') 'mean abs(mag) per spin = ', abs_mave
      WRITE(*,'(A28,F10.4)') 'mean squared mag. per spin = ', m2ave
      END

      SUBROUTINE save_config(N,L,T,spin)
      IMPLICIT NONE
      DOUBLE PRECISION T
      INTEGER N,L,x,y,spin(32,32)
      CHARACTER*10 filename
      WRITE(*,*)'name of file for last configuration = '
      READ(*,*) filename
      OPEN(2,FILE=filename,STATUS='NEW')
      WRITE(2,*) T
      DO y = 1,L
         DO x = 1,L
            WRITE(2,*) spin(x,y)
         END DO
      END DO
      CLOSE(2)
      END

      SUBROUTINE read_config(N,L,T,spin)
      IMPLICIT NONE
      DOUBLE PRECISION T
      INTEGER N,L,x,y,spin(32,32)
      CHARACTER*10 filename
      WRITE(*,*)'filename ? '
      READ(*,*) filename
      OPEN(1,FILE=filename,STATUS='OLD')
      READ(1,*)T
      DO y = 1,L
         DO x = 1,L
            READ(1,*) spin(x,y)
         END DO
      END DO
      CLOSE(2)
      END

      SUBROUTINE correl(Ce,Cm,E,M,esave,msave,pass,tsave)
*     accumulate data for time correlation functions
*     save last tsave values of M and E
      IMPLICIT NONE
      INTEGER tdiff,tsave,pass,index0,index
      DOUBLE PRECISION E,M,Ce(0:20),Cm(0:20),esave(100),msave(100)
*     index0 = array index for earliest saved time
      index0 = mod(pass-1,tsave) + 1
      IF (pass .GT. tsave) THEN
*        compute Ce and Cm after tsave values are saved
         index = index0
         DO tdiff = tsave,1,-1
            Ce(tdiff) = Ce(tdiff) + E*esave(index)
            Cm(tdiff) = Cm(tdiff) + M*msave(index)
            index = index + 1
            IF (index .GT. tsave) index = 1
         END DO
      END IF
*     save latest value in array position of earliest value
      esave(index0) = E
      msave(index0) = M
      END

      SUBROUTINE c_output(Ce,Cm,cum,mcs,tsave)
*     compute time correlation functions
      IMPLICIT NONE
      INTEGER tdiff,tsave,mcs
      DOUBLE PRECISION Ce(0:20),Cm(0:20),cum(10)
      DOUBLE PRECISION ebar,e2bar,mbar,m2bar,norm
      ebar = cum(1)/mcs
      e2bar = cum(2)/mcs
      Ce(0) = e2bar - ebar*ebar
      mbar = cum(3)/mcs
      m2bar = cum(4)/mcs
      Cm(0) = m2bar - mbar*mbar
      norm = 1.0D0/(mcs-tsave)
      WRITE(*,*)
      WRITE(*,*) '    t  ','  Ce(t) ','   Cm(t)'
      WRITE(*,*)
      DO tdiff = 1,tsave
*        define correlation function so that C(0) = 1
*        and C(infinity) = 0
         Ce(tdiff) = (Ce(tdiff)*norm - ebar*ebar)/Ce(0)
         Cm(tdiff) = (Cm(tdiff)*norm - mbar*mbar)/Cm(0)
         WRITE(*,15) tdiff,Ce(tdiff), Cm(tdiff)
15       FORMAT(I6,2(F10.5))
      END DO
      END

      REAL FUNCTION rnd(seed)
*     linear congruential random number generator with shuffling
*     based on ran1 in "Numerical Recipes" second edition
      IMPLICIT NONE
      INTEGER a,m,q,p,n,ndiv,j,k,seed
      REAL rm,rmax
*     m = 2**31 - 1  and  m = a*q + p
      PARAMETER (a = 16807, m = 2147483647, rm = 1.0/m)
      PARAMETER (q = 127773, p = 2836, n = 32, ndiv = 1 + (m-1)/n)
      PARAMETER (rmax = 1.0 - 1.2e-7)
      INTEGER r(n),r0,r1
      SAVE r,r0,r1
      DATA r/n*0/
*     initialize table of random numbers
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
*     beginning when not initializing
*     compute r1 = mod(a*r1,m) without overflows
      k = r1/q
      r1 = a*(r1 - k*q) - p*k
      IF (r1 .LT. 0) r1 = r1 + m
      j = 1 + r0/ndiv
      r0 = r(j)
      r(j) = r1
      rnd = min(rm*r0,rmax)
      END
