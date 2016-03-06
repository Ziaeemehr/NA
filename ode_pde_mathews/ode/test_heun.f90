program test_heun
implicit none
integer :: M,i
real(8) :: k1,k2,h,a,b,ya,f
real(8) :: t, y

ya= 1.d0
a = 0.d0
b = 3.d0
! Change M to modify the "h"
M = 24 
!# (M,h) = (3,1) , (6,0.5), (12,0.250),(24,0.125)
h = (b-a)/real(M)
y = ya
t = a
write (*,'(1x,A2,5x,A5)') "i","y"
write(*,100) t,y
do i = 1,M
    k1 = f(t,y)
    k2 = f(t + h, y + h * k1)
    y = y + 0.5 * h * ( k1 + k2 )
    write(*,100) t + h, y
    t = t+h
enddo
100 format (1x,f5.3,f10.6)
end program

function f(t,y)
    implicit none
    real(8),intent(in) :: t
    real(8),intent(in):: y
    real(8) :: f
    f =0.5*(t-y)
end function

