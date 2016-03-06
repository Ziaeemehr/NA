program test
implicit none
real(8) :: f,a,b,ya
integer :: M
interface
function f(t,y)
    real(8),intent(in):: t,y
end function f
end interface
ya= 0.d0
a = 0.d0
b = 3.d0
M = 30
write(*,*) "problem 9 page 473 Mathews Numerical methods using matlab 3rd Ed."
write(*,*) "RK4 Method"

call rk4(f,a,b,ya,M)
end program

function f(t,y)
implicit none
real(8),intent(in) :: t
real(8),intent(in):: y
real(8) :: f
    f = 1./exp(0.5 * t*t) 
end function

! OUTPUT 
!
!RK4 Method
! 0.000 0.0000000
! 0.500 0.6914625
! 1.000 0.8413448
! 1.500 0.9331928
! 2.000 0.9772499
! 2.500 0.9937903
! 3.000 0.9986501
