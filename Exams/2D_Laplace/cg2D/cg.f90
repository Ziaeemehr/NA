program lap_cg
! laplace equation with conjugate gradient algorithm 
! not saving the whole matrix
! to plot the result:
! echo 'splot "r.txt" nonuniform matrix w l lt 7' | gnuplot --persist
implicit none
integer :: N ,M
integer :: nv,k,j,i,iter
real(8) :: h,hsq,dot,gamma,oldrdr,alpha,a,tol
real(8),allocatable :: x(:),b(:), p(:), r(:), v(:)
real(8) :: xc, yc,pi
! -----------------------------------------------
a = 100.d0
!b = 1.d0 
M = 1000
tol = 1.d-5
pi = 4.d0 *atan(1.d0)
h = a/dble(M);
hsq = h*h;
N = M/100
nv = (N-1)*(M-1)
print *, M,N

allocate(x(0:nv-1),b(0:nv-1), p(0:nv-1), r(0:nv-1), v(0:nv-1))

! -----------------------------------------------
! build vector b
k = 0
do i = 1,M-1
    xc = i*h
    do j = 1,N-1
       yc = j*h       
       b(k) = 0.d0
!       if (i==1) then
!            !left edge -- include -u(0,y)/h^2 in RHS
!            b(k) = b(k) - yc*yc/hsq
!       elseif( i == (N-1) ) then
!            !right edge -- include -u(1,y)/h^2 in RHS
!           b(k) = b(k) - (yc*yc-1.d0)/hsq
!       endif
       if( j==1 ) then
           ! bottom edge -- include -u(x,0)/h^2 in RHS
!           b(k) = b(k) - xc*xc/hsq
       elseif( j==(N-1) ) then
           ! top edge -- include -u(x,1)/h^2 in RHS
            b(k) = b(k) - (sin(pi*xc/a))/hsq       
       endif
       k = k+1
    enddo
enddo
! -----------------------------------------------
! zero x to start
x(0:nv-1) = 0.d0
r(0:nv-1) = b(0:nv-1)
p(0:nv-1) = b(0:nv-1)

iter = 0

do  while (1)
    iter = iter +1
    dot = 0.d0
    do i = 0,nv-1
        dot = dot + p(i) * p(i)
    enddo
    dot = sqrt(dot)
    write(*,'(I6,es18.9)') iter, dot
    if( dot < tol ) exit
    
    ! v = A*p
    call Ax( v, p, nv, N, hsq )

    oldrdr = 0.d0
    dot = 0.d0
    do i=0,nv-1
        oldrdr = oldrdr + r(i)*r(i)
        dot = dot + p(i)*v(i)
    enddo
     
    alpha = oldrdr/dot
    
    do i = 0,nv-1 
        x(i) = x(i) + alpha * p(i)         
    enddo
            
    dot = 0;    
    do i = 0,nv-1
       r(i) =r(i) - alpha*v(i)
       dot = dot + r(i)*r(i) 
    enddo
    
    gamma = - dot/oldrdr

    do i = 0,nv-1 
        p(i) = r(i) - gamma*p(i)
    enddo
enddo
open (unit=10, file = 'r.txt')
write(10,120) nv,(i*h,i=1,N-1)
!write(10, '(9f15.7)') x
do i=0,M-2
    write(10,'(*(f15.7))') (i+1)*h,x(0+i*(N-1):N-2+i*(N-1))
    !write(10,'(*(f15.7))') (i+1)*h,x(0+i*(M-1):M-2+i*(M-1))
enddo
    
    

!open (unit=11, file = 'x.txt')
!do i = 1,N-1
!write (11,'(f15.7)') i*h
!enddo


120 format (I15,*(f15.7))
deallocate(x,b,p,r,v)
  contains
! -----------------------------------------------
  ! compute p = Ax without storing A
subroutine Ax(p,x,nv,N,hsq)
integer :: nv,N,i,j
real(8) :: p(0:nv-1), x(0:nv-1)
real(8) :: hsq

do k = 0,nv-1
    j = k/(N-1)+1
    i = mod(k,(N-1)) + 1
    p(k) = -4.d0 * x(k)/hsq
    if( i/=1   ) p(k) = p(k) + x(k-1)/hsq;
    if( i/=N-1 ) p(k) = p(k) + x(k+1)/hsq;
    if( j/=1   ) p(k) = p(k) + x(k-(N-1))/hsq
    if( j/=N-1 ) p(k) = p(k) + x(k+(N-1))/hsq
enddo
end subroutine

end program
