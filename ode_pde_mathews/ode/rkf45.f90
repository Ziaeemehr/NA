subroutine rkf45(f,a,b,ya,M)
   
implicit none
integer :: i, M
real(8) :: a2,a3,a4,a5,a6
real(8) :: b2,b3,b4,b5,b6
real(8) :: c3,c4,c5,c6
real(8) :: d4,d5,d6,e5,e6
real(8) :: r1,r3,r4,r5,r6
real(8) :: f6, h
real(8) :: n1,n3,n4,n5
real(8) :: m1,m3,m4,m5,m6
real(8) :: a,b,k1,k2,k3,k4,k5,k6,ya
real(8) :: t,y,f
real(8) :: delta,R,eps,y1,y2

interface 
function f(t,y)
    real(8),intent(in):: t,y
end function
end interface


h = (b-a)/real(M)
y = ya
t = a
eps=1.d-7
a2=1./4; b2=a2; a3=3./8; b3=3./32; c3=9./32; a4=12./13;
b4=1932./2197; c4=-7200./2197; d4=7296./2197; a5=1;
b5=439./216; c5=-8; d5=3680./513; e5=-845./4104; a6=0.5;
b6=-8./27; c6=2.; d6=-3544./2565; e6=1859./4104; f6=-11./40; 
!r1=1./360; r3=-128./4275; r4=-2197./75240; r5=1./50;r6=2./55;
n1=25./216; n3=1408./2565; n4=2197./4104; n5=-0.2; 
m1=16./135; m3=6656./12825; m4=28561./56430; m5=-9./50; m6=2./55;
i = 0
write (*,'(A2,A7,A7,A8)') "i","t","y","h"
write(*,100) i,t,y,h

do while (t<b) 
    h = min(h,b-t)
    k1 = h * f(t,y)
    k2 = h * f(t + a2 * h, y + b2 * k1 )
    k3 = h * f(t + a3 * h, y + b3 * k1 + c3 * k2 )
    k4 = h * f(t + a4 * h, y + b4 * k1 + c4 * k2 + d4 * k3)
    k5 = h * f(t + a5 * h, y + b5 * k1 + c5 * k2 + d5 * k3 + e5 * k4)
    k6 = h * f(t + a6 * h, y + b6 * k1 + c6 * k2 + d6 * k3 + e6 * k4 + f6 * k5)
    y1 = y + n1 * k1 + n3 * k3 + n4 * k4 + n5 * k5
    y2 = y + m1 * k1 + m3 * k3 + m4 * k4 + m5 * k5 + m6 * k6
    R = abs(y1-y2)/h
    delta = 0.84 * (eps/R)**(1./4)
    if (R <= eps) then
        t = t + h
        y = y1
        i = i+1
        write(*, 100) i,t,y,h
        h = delta * h
    else
        h = delta * h
    endif
enddo
100 format (I2,f7.4,f10.7,f7.4)
end subroutine
