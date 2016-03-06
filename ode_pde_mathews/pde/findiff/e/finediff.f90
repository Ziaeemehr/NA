
! Finite difference Solution for the Wave Equation.
! To approximate the solution of 
!          utt(x,t) = c^2 uxx(x,t)
! over R={(x,t): 0<= x <= a,  0<= t <= b}
! with u(0,t)=0, u(a,t)= 0 and
! u(x,0)  = f(x)
! ut(x,0) = g(x) 
! 

program test_finediff
implicit none
integer :: i,j
integer, parameter :: n = 10
integer, parameter :: m = 10
real(8) :: a,b,s,c
real(8), dimension(n,m) :: u

interface 
   function f(x)
     real(8) :: x
   end function f
   function g(x)
     real(8) :: x
   end function g
end interface
 c = 2.d0
 a = 1.d0
 b = 5.d-1
 u = 0.d0
 call finediff(f,g,a,b,c,n,m,u)

 do i = 1,n
     do j=1,m
         write(*,100) (u(i,j))
     enddo
 enddo
 100 format (10f15.6)
end program
! ------------------------------------40
function f(x)
implicit none 
real(8) :: x
real(8) :: pi
real(8) :: f
pi = 4.d0 * atan(1.d0)
f = sin(pi * x) + sin(2.d0 * pi * x)
end function
! ------------------------------------40
function g(x)
implicit none
real(8) :: x
real(8) :: g
g = 0.d0
end function
! ------------------------------------40
subroutine finediff(f,g,a,b,c,n,m,u)
integer :: i,j
integer :: n 
integer :: m 
real(8) :: a,b,h,k,x,s,c
real(8), dimension(n,m) :: u

interface 

function f(x)
real(8) :: x
end function f

function g(x)
real(8) :: x
end function g

end interface

 h = a/(n-1)
 k = b/(m-1)
 c = 2.d0
 r = c*k/h
 r2 = r * r
 s = 1.d0 - r2
 
do i = 2, n-1
    x = h * (i-1)
    u(i,1) = f(x)
    u(i,2) = s * f(x) + k * g(x) + 5.d-1 * r2 * (f(x+h)+f(x-h))
enddo

do j = 3,m
    do i = 2,n-1
        u(i,j) = 2.d0 * s * u(i,j-1) + r2 * &
        (u(i-1,j-1) + u(i+1,j-1)) - u(i,j-2)
    enddo
enddo

end subroutine