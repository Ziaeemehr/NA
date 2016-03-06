!  schro - Program to solve the Schrodinger equation
!  for a free particle using the Crank-Nicolson scheme

program schro
 implicit none
 integer :: info,i,j,nStep,iStep
 integer,parameter :: N = 200
 double precision,parameter :: h_bar=1.d0, tau = 1.d0, mass = 1.d0
 double precision :: L,h,pi,x0,velocity,k0,sigma0      
 double precision :: Norm_psi,expFactor,probability(N),x(N)
 double complex :: Psi(N),PsiInit(N),NewPsi(N),pd(N),pa(N)
 double complex :: dl(N),dd(N),du(N),v(N),u(N),xx(N),y(N),q(N),temp(N)
 double complex :: iImag

 !* Initialize parameters (grid spacing, time step, etc.)      
 L = 100.d0
 h = L/dble(N-1)

 do i = 1,N
   x(i) = h*(i-1) - L/2   ! Coordinates  of grid points
 enddo

 !* Set up the Tridiagonal matrix plus 2 extra elements
 call MatrixQ(dl,dd,du,u,v,N,tau,h)
 
!* Initialize the wavefunction
 pi = 4.d0*atan(1.d0)
 x0 = 0.d0                                            ! Location of the center of the wavepacket
 velocity = 0.5                                       ! Average velocity of the packet
 k0 = mass*velocity/h_bar                             ! Average wavenumber
 sigma0 = L/10.d0                                     ! Standard deviation of the wavefunction
 Norm_psi = 1/(sqrt(sigma0*sqrt(pi)))                 ! Normalization
 iImag = ( 0.0, 1.0 )    ! = sqrt(-1)
 do i=1,N
   expFactor = exp(-(x(i)-x0)*(x(i)-x0)/(2*sigma0*sigma0))
   Psi(i) = Norm_psi * cos(k0*x(i)) * expFactor &
      + iImag* Norm_psi * sin(k0*x(i)) * expFactor
 enddo
 
!  nStep = int(L/(velocity*tau))  ! Particle should circle system
!  print*, nStep
 nStep =40
 
 do i=1,N                       ! Record initial condition
   probability(i) = abs(Psi(i)*Psi(i))
 enddo
 
 ! printing the psi and probability 
 open(unit=7,file='realpsi.txt')
 open(unit=8,file='imgpsi.txt')
 open(unit=9,file='prob.txt')
 write(7,110) N,(x(i),i=1,N)
 write(7,110) 0,(real(Psi(i)),i=1,N)         ! at t=0
 write(8,110) N,(x(i),i=1,N)
 write(8,110) 0,(AIMAG(Psi(i)),i=1,N)        ! at t=0
 write(9,110) N,(x(i),i=1,N)
 write(9,110) 0,(probability(i),i=1,N)       ! at t=0
 110 format (I15,*(es15.7))
!  120 format (*(f15.7))
   
 !* Loop over desired number of steps (wave circles system once)
 do iStep = 1,nStep 
   temp = Psi
   call zgtsv( N, 1 , dl , dd , du , temp , N , info )
   y = temp
   call MatrixQ(dl,dd,du,u,v,N,tau,h)   
   call zgtsv( N, 1 , dl , dd , du , u , N , info )
   q = u

   xx = y - dot_product(v,y)/(1.d0+dot_product(v,q))*q 
   NewPsi = xx - Psi
   probability =   abs(NewPsi * NewPsi)
   
   write(7,110)  iStep, (real(NewPsi(i)), i=1,N)
   write(8,110)  iStep, (AIMAG(NewPsi(i)),i=1,N)
   if (mod(iStep,20)==0) then
     write(9,110)  iStep, (probability(i),  i=1,N)
   endif
   
   Psi = NewPsi   
   call MatrixQ(dl,dd,du,u,v,N,tau,h)
 enddo
 
 contains
 
 subroutine MatrixQ(dl,dd,du,u,v,N,tau,h)
 integer :: N
 double precision :: coeff,h_bar,tau,masss,h
 double complex :: iImag,pa,pd
 double complex :: dl(N),dd(N),du(N),u(N),v(N)
 !* Set up the Tridiagonal matrix plus 2 extra elements
 
 h_bar =1.d0
 coeff = -h_bar**2/(2*mass*h**2)
 iImag = ( 0.0, 1.0 )    ! = sqrt(-1)
 pd = 0.5 * (1.d0 - 2.d0 * coeff * iImag*0.5*tau/h_bar)
 pa = 0.5 * coeff * iImag*0.5*tau/h_bar
 dl = pa
 dd = pd
 du = pa
 dd(1) = 2.d0*pd
 dd(N) = pd+pa*pa/pd
 
 ! Define tow Sherman-Morison vectors
 u = 0.d0
 u(1) = - pd
 u(n) = pa
 v = 0.d0
 v(1) = 1.d0
 v(n) = -pa/pd
   
 end subroutine
 
end program  

!         