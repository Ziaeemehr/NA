subroutine rk2(f,a,b,ya,M)
!
! (Runge-Kutta Method of Order 2)
! To approximate the solution of the inital value problem y' = f(t,y)
! with  y(a) = y0 over [a,b] by using formula:
! y k+1 = y k + h*f(t + h/2, yk + h/2 * f(t,yk))

implicit none
integer :: i,M
real(8) :: k1,k2,y,t,ya,h,a,b
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
    y  = y +  k2
    t  = t + h 
    write (*,100) t, y
enddo
100 format (1x,f5.3,f10.6)
end subroutine
