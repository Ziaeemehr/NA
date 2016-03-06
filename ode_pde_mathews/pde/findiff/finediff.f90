
! Finite difference Solution for the Wave Equation.
! To approximate the solution of 
!          utt(x,t) = c^2 uxx(x,t)
! over R={(x,t): 0<= x <= a,  0<= t <= b}
! with u(0,t)=0, u(a,t)= 0 and
! u(x,0)  = f(x)
! ut(x,0) = g(x) 
! plot the result
! echo 'splot "result.txt" nonuniform matrix w l lt 7' | gnuplot --persist
program test_finediff
implicit none
integer :: i,j
integer, parameter :: n = 11
integer, parameter :: m = 11
real(8) :: a,b,h,k,x,c,r,f
real(8), dimension(n,m) :: u

 c = 2.d0
 a = 1.d0
 b = 5.d-1
 u = 0.d0
 
 h = a/(n-1)
 k = b/(m-1)
 c = 2.d0
!  r = c*k/h
 r = 1.d0
 
do i = 2, n-1
    x = (i-1)* h
    u(i,1) = f(x)
    u(i,2) = 0.5 * (f(x+h)+f(x-h))
enddo

do j = 3,m
    do i = 2,n-1
        u(i,j) =  u(i+1,j-1) + u(i-1,j-1) - u(i,j-2)
    enddo
enddo
open(unit=7,file='result.txt')

write(7,110) n*m,((i-1)*h,i=1,n)
110 format (I12,*(f12.6))

do j = 1,m 
    write(7,120) (j-1)*k, (u(i,j),i=1,n)
enddo

120 format (*(f12.6))
end program
! -----------------------------------------------
function f(x)
implicit none 
real(8),intent(in) :: x
real(8) :: pi
real(8) :: f

pi = 4.d0 * atan(1.d0)
f = sin(pi * x) + sin(2.d0 * pi * x)

end function

