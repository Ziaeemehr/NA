program test_ddot
implicit none
!include 'mkl.fi'
!include 'blas.f90'
integer,parameter :: n = 2
real(kind=8) :: x(n), y(n)
integer :: incx, incy
real(kind=8) :: res,res2
real(kind=8) :: ddot

incx = 1
incy = 1
x = (/1.d0 ,2.d0/)
y = (/2.d0, 3.d0/)

res = ddot (n,x,incx,y,incy)
!res2 = dot_product(x,y)
write(*,*) res
end program
