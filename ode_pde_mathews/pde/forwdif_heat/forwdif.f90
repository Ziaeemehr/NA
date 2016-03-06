program forwdiff
! Forward-difference Method for Heat Equation.
! To approximate the solution of heat equation 
!    Ut(x,t) = c * Uxx(x,t)
! x in [0,a] , t in [0,b] 
! U(x,0) = f(x), U(0,t) = c1 , U(a,t) = c2
! n and m number of grid points over [0,a] and [0,b]

implicit none
integer,parameter :: nx = 6 , mt = 11
real(8) :: a,b,c,h,k,r,s,x!,fun
real(8) :: U(nx,mt)
integer :: i,j
!interface
!   function f(x)
!   real(8),intent(in) :: x
!   end function f
!end interface
 c = 1.d0
 b = 2.d-1 
 a = 1.d0
 h = a/real(nx-1)
 k = b/real(mt-1)
 r = c * c * k /(h * h)

s = 1.d0 - 2.d0 * r
U (1:nx,1:mt) = 0.d0

! Boundary condition
U(1,1:mt) = 0.d0
U(nx,1:mt) = 0.d0

!Generate the first row
do i = 2,nx-1
    x = (i-1) * h
    U(i,1) = fun(x)
enddo

! Generate remaining rows of U
do j = 2,mt
    do i = 2,nx-1
        U(i,j) = s * U(i,j-1) + r * (U(i-1,j-1) + U(i+1,j-1))
    enddo
enddo
write(*,120) ((i-1)*h,i=1,nx)
120 format (7x,' t',7x,'x = ',f6.2,5f12.2)
do j = 1,mt
    write(*,110) (j-1)*k, (U(i,j),i=1,nx)
enddo
110 format (f12.2,6f12.6)

 contains
function fun(x)
    implicit none
    real(8) :: x
    real(8) :: fun
    fun = 4.d0 * x - 4.d0 * x * x
end function fun
end program   
   