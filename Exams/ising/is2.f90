program ising2D
implicit none
integer, parameter :: seed = 123459
integer, parameter :: L = 32
integer :: N,i,j,pass,nequil
integer :: spin(L,L)
real(8) :: E,M,T,w(-8:8)
real(8) :: rnd,r,accept



T = 2.D0
N = L*L
E = 0.d0
M = 0.d0
nequil = L*L    !MC steps per spin for equilibrium

CALL initial(L,T,E,M,w,spin)

! print *, w
! stop
DO pass = 1,nequil
    CALL Metropolis(N,L,spin,E,M,w,accept)
END DO

! do i =1,L
!     write(*,'(*(I5))') (spin(i,j),j=1,L)
! enddo




end program

SUBROUTINE initial(L,T,E,M,w,spin)
IMPLICIT NONE
integer :: spin(L,L)
integer :: L,x,y,up,right,seed,dE
real :: rnd,dummy
real(8) :: T,E,M,summ,w(-8:8)

!        seed must not equal 0 to initialize rnd
dummy = rnd()

!        random initial configuration
DO y = 1,L
    DO x = 1,L       
       IF (rnd() .LT. 0.5) THEN
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
! print *, E,M
!        compute Boltzmann probability ratios
DO dE = -8,8,4
   w(dE) = exp(-dE/T)
END DO
!print *, w
!tsave = 10
END SUBROUTINE
! ---------------------------------------------------------


subroutine Metropolis(N,L,spin,E,M,w,accept)
!  one Monte Carlo step per spin
   integer :: ispin,x,y,dE,L
   integer :: spin(L,L)
   real :: rnd
   real(8) :: w(-8:8),E,M,accept
   do ispin = 1,N
!     random x and y coordinates for trial spin
      call random_number(rnd)
      x = int(L*rnd) + 1
      call random_number(rnd)
      y = int(L*rnd) + 1
      dE = DeltaE(x,y,L,spin)
      print *,dE      
      call random_number(rnd)
      if (rnd <= w(dE)) then
         spin(x,y) = -spin(x,y)
         accept = accept + 1
         M = M + 2*spin(x,y)
         E = E + dE
      end if
   end do
end subroutine metropolis

INTEGER FUNCTION DeltaE(x,y,L,spin)
! periodic boundary conditions
IMPLICIT NONE
INTEGER L,x,y,left,right,up,down,spin(L,L)
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

! ---------------------------------------------------------
FUNCTION rnd()
real(8) :: rnd,x
call random_number(x)
rnd = x
 END FUNCTION
! ---------------------------------------------------------
