program test_traprl
! composite Tramzoidal Rule
implicit none
real(8) :: a,b,s
integer :: M

a = 1.d0
b = 6.d0
M = 10
s = 0.d0
call traprl(f,a,b,M,s)
write(*,'(f12.8)') s

 contains
!-------------------------------------- 
function f(x)
real(8),intent(in) :: x
real(8) :: f
f = 2.d0 + sin(2.d0 * sqrt(x))
end function
!--------------------------------------
subroutine traprl(f,a,b,M,s)
! To approximate the integral f(x) over [a,b]
! using composite Tramzoidal Rule
! x(k) = a+k*h k=0,...M, x(0) = a , x(M)= b
! Numerical methods, Mathews page 363
implicit none
real(8) :: h,a,b,s,x,f
integer :: M,k

h = (b-a)/real(M)
s = 0.d0
do  k = 1, (M-1)
    x = a + h * k
    s = s + f(x)
enddo

s = h * (f(a)+f(b)) * 5.d-1 + h * s
end subroutine
!--------------------------------------
end program
