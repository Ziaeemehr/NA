!Garsia chapter 8
! relax - Program to solve the Laplace equation using
! Jacobi, Gauss-Seidel and SOR methods on a square grid
!f1,f2,f3,f4 : botom,top,lef and right faces respectively
program sor
integer,parameter :: MAXN=300
integer :: Nx,Ny,i,j,iterMax,iter,nIter
real(8) :: Lx,Ly,x,y,h,omega,omegaOpt,pi
real(8) :: coeff,phiTemp,ave
real(8) :: changeDesired, change, changeSum
real(8),allocatable :: phi(:,:)


pi = 4.d0*atan(1.d0)
Nx = 9
Ny = 9
Lx = 4.d0             ! System size (length)
Ly = 4.d0
h = Lx/dble(Nx-1)       ! Grid spacing
! print *, h

allocate(phi(Nx,Ny))


r = 5.d-1*(cos(pi/Nx)+cos(pi/Ny))
omega = 2.d0/(1.d0+sqrt(1-r*r))
print *, "omega", omega

! set initial guess
ave = (Ly *(f1(0.d0) + f2(0.d0)) + Lx * (f3(0.d0) + f4(0.d0)))/(2.d0*Lx + 2.d0*Ly)
phi(1:Nx,1:Ny) = ave

!* Set boundary conditions
do j = 1,Ny
    y = (j-1) * h
    phi(1,j)  = f3(y)
    phi(Nx,j) = f4(y)
enddo
do i = 1,Nx
    x = (i-1) * h
    phi(i,1)  = f1(x)
    phi(i,Ny) = f2(x)
enddo


phi(1,1) = (phi(1,2)   + phi(2,1))  * 5.d-1
phi(1,Ny) = (phi(1,Ny-1) + phi(2,Ny))  * 5.d-1
phi(Nx,1) = (phi(Nx-1,1) + phi(Nx,2))  * 5.d-1
phi(Nx,Ny) = (phi(Nx-1,Ny) + phi(Nx,Ny-1))* 5.d-1
!* Loop until desired fractional change per iteration is obtained
 iterMax = Nx*Ny             ! Set max to avoid excessively long runs
 changeDesired = 1e-4        ! Stop when the change is given fraction
write(*,*) 'Desired fractional change = ', changeDesired
do iter=1,iterMax
  changeSum = 0
  do i=2,(Nx-1)            ! Loop over interior points only
    do j=2,(Ny-1)
      phiTemp = 0.25*omega*(phi(i+1,j)+phi(i-1,j) &
      + phi(i,j-1)+phi(i,j+1))  +  (1-omega)*phi(i,j)
      changeSum = changeSum + abs(1-phi(i,j)/phiTemp)
      phi(i,j) = phiTemp
    enddo
  enddo
  
  !* Check if fractional change is small enough to halt the iteration
  change = changeSum/((Nx-2)*(Ny-2))
  write(*,*) iter,change
  if( change .lt. changeDesired ) then
    write(*,*) "Desired accuracy achieved after ", iter," iterations"
    write(*,*) "Breaking out of main loop"
    nIter = iter
    exit       ! Break out of the main loop
  endif
enddo

open(13,file='phi.txt')

write(13,110) Nx*Ny,((i-1)*h,i=1,Nx)
do j = 1,Ny
  write(13,120) (j-1)*h,(phi(i,j),i=1,Nx)
enddo

110 format(I12,*(f12.6))
120 format(*(f12.6))
deallocate(phi)
 contains
 function f1(x)
     real(8),intent(in) :: x
     real(8) :: f1
     f1 = 20.d0
 end function f1
 function f2(x)
     real(8),intent(in) :: x
     real(8) :: f2
     f2 = 180.d0
 end function f2
 
 function f3(x)
     real(8),intent(in) :: x
     real(8) :: f3
     f3 = 80.d0
 end function f3
 
 function f4(x)
     real(8),intent(in) :: x
     real(8) :: f4
     f4 = 0.d0
 end function f4
end program
