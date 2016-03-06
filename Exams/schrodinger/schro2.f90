!  schro - Program to solve the Schrodinger equation
!  for a free particle using the Crank-Nicolson scheme

program schro
 implicit none
 integer :: info,i,j,nStep,iStep
 integer,parameter :: N = 500
 double precision,parameter :: h_bar=1.d0, tau = 1.d0, mass = 1.d0
 double precision :: L,h,pi,x0,velocity,k0,sigma0      
 double precision :: Norm_psi,expFactor,probability(N),x(N)
 double complex :: Psi(N),PsiInit(N),NewPsi(N),pd(N),pa(N),temp(N)
 double complex :: QQ(N,N)
 double complex :: iImag
 integer :: ipiv(N)

 !* Initialize parameters (grid spacing, time step, etc.)      
 L = 100.d0
 h = L/dble(N-1)
 print *, h
 do i = 1,N
   x(i) = h*(i-1) - L/2   ! Coordinates  of grid points
 enddo
 
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
 
 nStep = int(L/(velocity*tau))  ! Particle should circle system
 nStep =10
!  do i=1,N                       ! Record initial condition
!    probability(i) = abs(Psi(i)*Psi(i))
!  enddo
 
 ! printing the psi and probability 
 open(unit=7,file='realpsi.txt')
!  open(unit=8,file='prob.txt')
 write(7,110) N,(x(i),i=1,N)
 write(7,110) 0,(real(Psi(i)),i=1,N)        ! at t=0
 110 format (I15,*(es15.7))

 call MatrixQQ(QQ,N,tau,h)

 
 do iStep = 1,nStep 
   temp = Psi
   call zgesv ( N , 1 , QQ , N , ipiv , temp , N , info )
   NewPsi =  temp - Psi
   write(7,110)  iStep, (real(NewPsi(i)),i=1,N)
   Psi = NewPsi   
   call MatrixQQ(QQ,N,tau,h)
 enddo
 
 contains
  
 subroutine MatrixQQ(QQ,N,tau,h)
 integer :: N
 double precision :: coeff,h_bar,tau,masss,h
 double complex :: iImag,pa,pd
 double complex :: dl(N),dd(N),du(N),u(N),v(N),QQ(N,N)
 !* Set up the Tridiagonal matrix plus 2 extra elements
 
 h_bar =1.d0
 coeff = -h_bar**2/(2*mass*h**2)
 iImag = ( 0.0, 1.0 )    ! = sqrt(-1)
 pd = 0.5 * (1.d0 - 2.d0 * coeff * iImag*0.5*tau/h_bar)
 pa = 0.5 * coeff * iImag*0.5*tau/h_bar
 dl = pa
 dd = pd
 du = pa
 do i=1,N
   do j = 1,N
   if (i==j)     QQ(i,j) = pd
   if (i==j+1)   QQ(i,j) = pa
   if (i+1==j)   QQ(i,j) = pa
   enddo
 enddo
 QQ(1,N) = pa
 QQ(N,1) = pa
 
 end subroutine
 
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