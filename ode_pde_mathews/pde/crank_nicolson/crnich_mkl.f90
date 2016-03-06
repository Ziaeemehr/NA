program crnich
!
! Crank-Nicolson Method for the Heat Equation.
! ! To approximate the solution of heat equation 
!    Ut(x,t) = c * Uxx(x,t)
! x over [0,a] , t over [0,b] 
! U(x,0) = f(x), U(0,t) = c1 , U(a,t) = c2
! n and m number of grid points over [0,a] and [0,b]

implicit none
integer,parameter :: nx = 11 , mt = 11
real(8) :: a,b,c,h,k,r,s1,s2,x,c1,c2
real(8) :: U(nx,mt),ux(nx),vd(nx),vb(nx),vc(nx),va(nx)
real(8) :: mat(nx,nx),temp(nx,nx)
integer :: i,j,info
integer, dimension(nx) :: ipiv
 c1 = 0.d0
 c2 = 0.d0
 c  = 1.d0
 b  = 1.d-1    ! t over [0,b]
 a  = 1.d0     ! x over [0,a]
 h  = a/real(nx-1)
 k  = b/real(mt-1)
 r  = 1
 s1 = 4.
 s2 = 0.

U (1:nx,1:mt) = 0.d0

! Boundary condition
U(1,1:mt)  = 0.d0
U(nx,1:mt) = 0.d0

!Generate the first row
do i = 2,nx-1
    x = (i-1) * h
    U(i,1) = f(x)
enddo
! Form diagonal and off-diagonal elements of A and
! the constant vector B and solve tridiagonal system Ax=B
!vd(1:nx) = s1
!vd(1) = 1.d0
!vd(nx) = 1.d0
!va(1:nx) = -1.d0
!va(nx-1) = 0.d0
!vc(1:nx) = -1.d0
!vc(1) = 0.d0
vb(1:nx) = 0.d0
vb(1) = 0.
vb(nx)= 0.
mat = 0.d0
do i =2,nx
    do j = 2,nx
    if (i==j)      mat(i,j) = 4.d0
    if (i==(j+1))  mat(i,j) = -1.d0
    if (j==(i+1))  mat(i,j) = -1.d0
    enddo
enddo
mat(nx,nx) = 1.d0;   mat(1,2) = 0.d0
mat(nx,nx-1) = 0.d0; mat(1,1) = 1.d0
temp = mat
do j = 2,mt    
    do i = 2,nx-1
        vb(i) = U(i-1,j-1) + U(i+1,j-1)+ s2 * U(i,j-1)
    enddo
    !write (*,"(20f10.4)") vb
    call dgesv(nx, 1, temp, nx, ipiv, vb, nx, info)
    U(1:nx,j) = vb
    temp = mat
enddo

! Printing the results
write(*,120) ((i-1)*h,i=2,nx-1)
120 format (7x,' t',7x,'x = ',f5.2,9f11.2)
do j = 1,mt
    write(*,110) (j-1)*k, (U(i,j),i=2,nx-1)
enddo
110 format (f12.2,999f12.6)

! Exact result for U(x,t)
! f = sin(pi*x)* 1.d0/exp(pi*pi*t)+sin(3*pi*x)*1.d0/exp(-9.d0*pi*pi*t)
 
 contains
 
function f(x)
    implicit none
    real(8) :: x
    real(8) :: f
    real(8) :: pi
    pi = 4.d0 * atan(1.d0)
    f = sin (pi*x) + sin(3*pi*x)    
end function f

end program
