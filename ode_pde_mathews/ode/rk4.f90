subroutine rk4(f,a,b,ya,M)
!
! To approximate the solution of the inital value problem y' = f(t,y)
! with  y(a) = y0 over [a,b] by using formula:
! y k+1 = y k + h/6(k1+ 2*k2 + 2*k3 + k4)

implicit none
integer :: i,M
real(8) :: k1,k2,k3,k4,y,t,ya,h,a,b
real(8) :: f
interface
function f(t,y)
    real(8),intent(in):: t,y
end function f
end interface
h = (b-a)/real(M)
t = a
y = ya
write(*,100) t,y
do i = 1,M
    k1 = h * f(t,y)
    k2 = h * f(t + 0.5 * h, y + 0.5 * k1)
    k3 = h * f(t + 0.5 * h, y + 0.5 * k2)
    k4 = h * f(t + h, y + k3)
    y  = y + (k1 + 2 * k2 + 2 * k3 + k4)/6.
    t  = t + h 
    write (*,100) t, y
enddo
100 format (1x,f5.3,f12.9)
end subroutine
