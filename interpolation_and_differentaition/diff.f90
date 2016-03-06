program test_diff
implicit none
integer,parameter :: n = 10
real(8) :: x(n) ,y(n),dy(n-1),ddy(n-2),dyc(n-2)
integer :: I

x = (/.1, .2, .3, .4, .5, .6, .7,.8,.9,1./)
y = (/ 0.0550090,   0.9315435,   0.6957469,   0.7394335,   0.4756588, &
     0.0052795,   0.2315714,   0.0206270,   0.9565449,0.5379002 /)

call diff1(x,y,dy,n)
write(*,110) dy
! ---------------------------

call diff2(x,y,ddy,n)
write(*,110) ddy
print *, " Central difference o(h**2) "
! ---------------------------
call diffcenter(x,y,dyc,n)
write(*,110) dy


110  format (90f16.9)
contains
! ---------------------------
subroutine diffcenter(x,y,dy,n)
implicit none
real(8) :: x(n),y(n)
real(8) :: dy(n-2)
integer :: n
real(8) :: h

h = x(2)-x(1)   
dy = (y(3:n)-y(1:n-2))/(2.d0 * h)

end subroutine
! ---------------------------
subroutine diff1(x,y,dy,n)
! forward difference first order
implicit none
real(8),intent(in) :: x(n),y(n)
real(8),intent(out):: dy(n-1)
integer :: n

dy = (y(2:n)-y(1:n-1))/(x(2:n)-x(1:n-1))

end subroutine
! ---------------------------
subroutine diff2(x,y,ddy,n)
implicit none
! central difference second order
real(8),intent(in) :: x(n),y(n)
real(8),intent(out):: ddy(n-2)
real(8) :: h
integer :: n

h= x(2)-x(1)
ddy= (y(3:n) - 2.d0 * y(2:n-1) + y(1:n-2))/(h*h)

end subroutine
! ---------------------------
end program



