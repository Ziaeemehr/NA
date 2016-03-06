! plot the result with
! echo 'splot "result.txt" nonuniform matrix w l lt 7' | gnuplot --persist
program cg_2d
  implicit none
  real(8) :: dx,dy,g_norm,pi,alph,tol,x
  real(8) :: gg,gn_new,S
  integer,parameter :: nx=100, ny=10
  real(8) :: g(nx,ny),g_new(nx,ny),u(nx,ny),h(nx,ny),Ah(nx,ny)
  real(8) :: EdotE(ny,ny)
  real(8),parameter :: a=100.d0 , b=1.d0
  integer :: i,j,iter
  
  tol = 1.d-6
  pi = 4.d0*atan(1.d0)
  dx = a/dble(nx-1)
  dy = b/dble(ny-1)
    
  u = 0.d0
  ! B.C
  ! -------------
  u(1,:)  = 0.d0
  u(nx,:) = 0.d0
  u(:,1)  = 0.d0
  do i = 1,nx
    x = (i-1)*dx
    u(i,ny) = sin(pi*x/a)
  enddo
  iter = 0
  call Ax(g,u,dx,dy,nx,ny)
  g_norm = innerprod(g,g,nx,ny)

  h = -g  
  do while (g_norm > tol)
    iter = iter + 1
    call Ax(Ah,h,dx,dy,nx,ny)    
    alph = g_norm/innerprod(h,Ah,nx,ny)
    u = u + alph * h
    g = g + alph * Ah
    gn_new = innerprod(g,g,nx,ny)    
    h = -g + (gn_new/g_norm) * h
    g_norm = gn_new
    print *, iter, g_norm
  enddo
    
  open (unit = 10, file = "result.txt")
  write(10,"(I10,*(f10.6))")  nx,((j-1)*dy,j=1,ny)

  do i=1,nx
    write(10,"(*(f10.6))") (i-1)*dx,(u(i,j),j=1,ny)
  enddo
 
  call E_field(EdotE,u,nx,ny,dx,dy)
  
  call traprl2D(EdotE,dx,dy,nx,ny,S)
  

  write (*,130) 5.d-1 * S
  130 format ("The Energy is : ", f14.7)


  
  contains
  subroutine traprl2D(f,dx,dy,nx,ny,S)
  ! composite Tramzoidal Rule integration
  implicit none
  real(8) :: dx,dy,a,b,c,d,x,y,f(nx,ny)
  real(8) :: S,s1,s2,s3,s4,s5
  integer :: nx,ny,i,j
  
  s = 0.d0;   s1 = 0.d0;  s2 = 0.d0;  s3 = 0.d0
  s4 = 0.d0;  s5 = 0.d0
  
  do i = 2, (nx-1)
      s1 = s1 + f(i,1)
      s2 = s2 + f(i,ny)
  enddo
  
  do j = 2,(ny-1)
      s3 = s3 + f(1,j)
      s4 = s4 + f(nx,j)
  enddo
  do j = 2,(ny-1)
      do i = 2,(nx-1)
          s5 = s5 + f(i,j)
      enddo
  enddo
  
  S = 1.d0/4.d0*dx*dy*(f(1,1) + f(nx,1) + f(1,ny) + f(nx,ny)&
      + 2.d0*(s1+s2+s3+s4)+4.d0*s5)
  end subroutine  
  ! -------------------------------------------------------
    subroutine E_field(EdotE,V,nx,ny,dx,dy)
  integer :: i,j,nx,ny
  real(8) :: Ex(nx,ny),Ey(nx,ny),V(nx,ny),EdotE(nx,ny),dx,dy
  
  do j = 1,ny
    do i = 1,nx
      if (i >(nx-2)) then
        Ex(i,j) = - (3.d0 * V(i,j) - 4.d0 * V(i-1,j) + V(i-2,j))/(2.d0 * dx)        
      else
        Ex(i,j) = - (-3.d0 * V(i,j) + 4.d0 * V(i+1,j) - V(i+2,j))/(2.d0 * dx)
      endif      
    enddo
  enddo
  
  do j = 1,ny
    do i = 1,nx
      if (j >(ny-2)) then
        Ey(i,j) = - (3.d0 * V(i,j) - 4.d0 * V(i,j-1) + V(i,j-2))/(2.d0 * dy)        
      else
        Ey(i,j) = - (-3.d0 * V(i,j) + 4.d0 * V(i,j+1) - V(i,j+2))/(2.d0 * dy)
      endif
    enddo
  enddo

  ! E.E = EdotE
  do j=1,ny
    do i = 1,nx
      EdotE(i,j) = (Ex(i,j)*Ex(i,j))+(Ey(i,j)*Ey(i,j))
    enddo 
  enddo
  end subroutine
  ! -------------------------------------------------------
  subroutine Ax(g,x,dx,dy,nx,ny)
  implicit none  
  real(8) :: x(nx,ny),g(nx,ny)
  real(8) :: dx,dy
  integer :: nx,ny,i,j
  
  g = 0.d0 
  do j = 2,ny-1
    do i = 2,nx-1
      g(i,j) =(x(i+1,j)-2.d0*x(i,j)+x(i-1,j))/(dx*dx)&
               +(x(i,j+1)-2.d0*x(i,j)+x(i,j-1))/(dy*dy)
    enddo
  enddo
  end subroutine
  ! -------------------------------------------------------
  function innerprod(p,q,nx,ny)
  implicit none
  real(8),dimension(nx,ny) :: p,q
  real(8):: innerprod
  integer :: i,j,nx,ny
  
  innerprod = 0.d0
  do j = 1,ny
    do i = 1,nx
      innerprod = innerprod + p(i,j)*q(i,j)
    enddo
  enddo
  end function
  ! -------------------------------------------------------

  end program