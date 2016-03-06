subroutine heun(f,a,b,ya,M)
!
! To approximate the solution of the initail value problem
! y'=f(t,y) with y(a)= y0 over [a,b] by computing 
! y k+1 = yk + h/2 * (f(tk,yk) + f(t k+1, yk + f(tk,yk)))
! for k = 0,..M-1
    
    
    
implicit none
integer :: M    
integer :: i
real(8) :: k1,k2,h,a,b,ya,f
real(8) :: t, y
interface 
function f(t,y)
    real(8),intent(in):: t,y
end function f
end interface
 
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

end subroutine
