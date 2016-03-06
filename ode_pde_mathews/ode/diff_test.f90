program test
implicit none
real(8) :: f,a,b,ya
integer :: M
interface
function f(t,y)
    real(8),intent(in):: t,y
end function f
end interface
ya= 1.d0
a = 0.d0
b = 3.d0
M = 24
!write(*,*) "Heun Method"
!call heun(f,a,b,ya,M)
!write(*,*) "RK2 Method"
!call rk2(f,a,b,ya,M)
write(*,*) "RK4 Method"
call rk4(f,a,b,ya,M)
!write(*,*) "RKF45"
!call rkf45(f,a,b,ya,M)
write(*,*) "ABM4"
call abm4(f,a,b,ya,M)

end program

function f(t,y)
implicit none
real(8),intent(in) :: t
real(8),intent(in):: y
real(8) :: f
    f =0.5*(t-y)
end function

!function f(t,y)
!implicit none
!real(8),intent(in) :: t
!real(8),intent(in):: y
!real(8) :: f,pi
!pi = 4.d0 * atan(1.d0)
!    f = 1./exp(0.5 * t*t) 
!end function

