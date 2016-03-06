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
integer :: i,j
 c1 = 0.d0
 c2 = 0.d0
 c  = 1.d0
 b  = 1.d-1    ! t over [0,b]
 a  = 1.d0     ! x over [0,a]
 h  = a/real(nx-1)
 k  = b/real(mt-1)
 r  = c*c *k/(h*h)
 s1 = 2 + 2/r
 s2 = 2/r - 2

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
vd(1:nx) = 4.
vd(1) = 1.d0
vd(nx) = 1.d0
va(1:nx) = -1.d0
!va(nx-1) = 0.d0
va(nx) = 0.d0       ! here it was the bug .va(nx-1) = 0.d0 should be va(nx) = 0.d0 since va(nx-1) is a row above the boundary condition
!vc(1:nx-1) = -1.d0
vc(1:nx) = -1.d0
vc(1) = 0.d0
vb(1:nx) = 0.d0
vb(1) = c1
vb(nx)= c2

do j = 2,mt    
    do i = 2,nx-1
        vb(i) = U(i-1,j-1) + U(i+1,j-1)+ s2 * U(i,j-1)
    enddo
    call tridag(va,vd,vc,vb,ux,nx)    
    U(1:nx,j) = ux
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
 
SUBROUTINE tridag(a,b,c,r,u,n)
INTEGER n,NMAX
REAL(8):: a(n),b(n),c(n),r(n),u(n)
PARAMETER (NMAX=500)
INTEGER j
REAL(8):: bet,gam(NMAX)

if(b(1).eq.0.)pause "tridag: rewrite equations"

bet=b(1)
u(1)=r(1)/bet
do j=2,n

        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.)pause "tridag failed"

        u(j)=(r(j)-a(j)*u(j-1))/bet
enddo 
        do  j=n-1,1,-1
                u(j)=u(j)-gam(j+1)*u(j+1)
        enddo 
return
END subroutine

function f(x)
    implicit none
    real(8) :: x
    real(8) :: f
    real(8) :: pi
    pi = 4.d0 * atan(1.d0)
    f = sin (pi*x) + sin(3*pi*x)    
end function f

end program
