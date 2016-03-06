subroutine abm4(f,a,b,ya,M)
implicit none
integer :: i,M
real(8) :: k1,k2,k3,k4,t,ya,h,a,b,w
real(8) :: f
real(8) :: y(M+1)

interface 
function f(t,w)
    real(8),intent(in) :: t,w
end function
end interface

t = a
y(1) = ya
w    = ya
h=(b-a)/real(M)
write (*,100) t, y(1)
do i = 1,3
    k1 = f(t,w)
    k2 = f( t + h/2, w + h/2 * k1) 
    k3 = f( t + h/2, w + h/2 * k2)
    k4 = f( t + h  , w + h   * k3)
    w = w + h/6. * (k1 + 2 * (k2 + k3) + k4 )
    t = t + h
    y(i+1) = w
    write (*,100) t, y(i+1)
enddo

do i = 4,M
    w = w + h/24. * (55. * f(t, y(i)) - 59. * f(t-h, y(i-1))+ 37. * & 
        f(t-2.*h, y(i-2)) - 9.* f(t-3.*h, y(i-3)))
    y(i+1) = y(i) + h/24. *(9. * f(t+h,w) + 19. * f(t,y(i)) - 5. *  &
        f(t-h,y(i-1)) + f(t-2.*h, y(i-2)))
        t = t + h
    write (*,100) t, y(i+1)
enddo
100 format (f7.4,f10.8)
end subroutine