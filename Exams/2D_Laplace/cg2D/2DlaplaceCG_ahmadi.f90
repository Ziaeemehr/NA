program conjugate_gradient_2D_laplace

! plot the result with
! echo 'splot "result.txt" nonuniform matrix w l lt 7' | gnuplot --persist
implicit none

real(8), allocatable :: u(:,:)
real(8), allocatable :: g(:,:)
real(8), allocatable :: h(:,:)
real(8), allocatable :: v(:,:)

real(8) :: dx,dy,alpha,beta,gnrm_new,gnorm,  pi 
integer , parameter :: nx = 100, ny = 10
integer :: it,i,j 
real(8) :: tolerance = 1.d-6  , a = 100.0, b = 1.0
real(8) :: inner_product 
 pi = 4.d0*atan(1.d0)


dx = a /(nx-1)
dy = b /(ny-1)

allocate(u(nx,ny),g(nx,ny),h(nx,ny),v(nx,ny))



u = 0.0                      !initial guess
!------------------B.C --------------------------
u(1,:) = 0.0
u(:,1) = 0.0
u(nx,:) = 0.0

do i = 1, nx
   u(i,ny) = sin(pi*(i-1)*dx/a)
enddo

!----------------------------------------------------
call residual(g,u,dx,dy,nx,ny)                    
h=g                                               
gnorm = inner_product(g,g,nx,ny)                 
it=0

do while (gnorm > tolerance )               
  it = it+1
  call matrixVectorProd(v,h,dx,dy,nx,ny)                
  alpha = gnorm/(inner_product(h,v,nx,ny))              
  u = u + alpha*h   
  g = g - alpha*v                                         
  gnrm_new = inner_product(g,g,nx,ny) 
  h = g + (gnrm_new/gnorm)*h                                 
  gnorm = gnrm_new
 write(*,*) it,gnorm
enddo

open (unit = 10, file = "result.txt")
write(10,"(I10,*(f10.6))")  nx,((j-1)*dy,j=1,ny)

do i=1,nx
write(10,"(*(f10.6))") (i-1)*dx,(u(i,j),j=1,ny)
enddo



end program conjugate_gradient_2D_laplace
!******************************************************************************
subroutine residual(g,u,dx,dy,nx,ny)      ! calculating gradient g
    implicit none
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    real(8), intent(in) :: dx,dy
    real(8), intent(in) :: u(nx,ny)
    real(8), intent(out):: g(nx,ny)

    integer :: i,j
 g = 0.0
   do j=2,ny-1
    do i=2,nx-1
       g(i,j) = - ( u(i+1,j) - 2.0*u(i,j) + u(i-1,j)) / (dx*dx)  - (u(i,j+1)-2.0*u(i,j) + u(i,j-1))/(dy*dy)
       
    enddo
    enddo


end subroutine residual

!***************************************************************************
function inner_product(a,b,n,m)

    implicit none
    integer, intent(in) :: n,m
    real(8), intent(in) :: a(n,m),b(n,m)

    integer i,j
    real(8) :: inner_product

    inner_product = 0.0

    do i=1,n
    do j=1,m
    inner_product = inner_product + a(i,j)*b(i,j)          !norm 
    enddo
    enddo

  
    
end function inner_product
!*********************************************************************************
subroutine matrixVectorProd(v,h,dx,dy,nx,ny)

    implicit none

    integer, intent(in) :: nx
    integer, intent(in) :: ny
    real(8), intent(in) :: dx,dy
    real(8), intent(in) :: h(nx,ny)

    real(8), intent(out) :: v(nx,ny)

    integer :: i,j
v =0.0
    do j=2,ny-1
      do i=2,nx-1
        v(i,j) = ( h(i+1,j) - 2.0*h(i,j) + h(i-1,j))/(dx*dx) &
        + ( h(i,j+1) - 2.0*h(i,j) + h(i,j-1))/(dy*dy)
      enddo
    enddo


end subroutine matrixVectorProd
!!!!!!!!!!!*************************************************************!!!!!!!!!!!!!!!!!!!!!!
