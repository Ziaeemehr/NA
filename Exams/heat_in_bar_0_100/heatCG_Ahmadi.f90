program heat_CG
implicit none

integer :: it, i, time, j
integer , parameter :: N = 1000 
real(8) :: T_old(N), g_new(N), h(N),  g(N), T_new(N)
real(8) :: dt , dx , k , x , aa ,summ, beta, gnrm ,alpha, dott

 dt = 2.d-1
 dx = 1.d0/N
 k = 0.2
 aa = (k*dt)/(dx**2)

!--------------  initial value  ----------------
 T_old = 0.0        
 T_old(1) = 100.0   
 T_old(N) = 0.0      

!-----------------------------------------

                   do time = 1, 100 

!--------------  calculate g  ----------------------
    T_new = 0.0                  ! initial guess for solution
   T_new(1) = 100.0
  T_new(N) = 0.0

    g(1) =  T_new(1) - T_old(1)
          do i = 2, N-1
           g(i) = -aa* T_new(i-1) +(2*aa + 1) * T_new(i) - aa* T_new(i+1) - T_old(i)
         enddo
    g(N) = T_new(N) - T_old(N) 
!--------------------------------------------------------
 ! write(*,*) g
!stop
   do  it = 1 , 2000
     if (it == 1)  h=-g

!----------------------calculating alpha -----------
 summ = 0.0
 summ = summ +(-aa)*(h(1)*h(2) + h(N-1)*h(N))
  do i = 2, N-2
      summ = summ + (-2*aa)* h(i)*h(i+1)
  enddo

 dott = 0.0
 dott = dott + h(1)**2 + h(N)**2
  do i= 2, N-1
         dott = dott +(2*aa+1)*(h(i)*h(i))
  enddo
  
  alpha = dot_product(g,g)/(dott+summ)      
      !  write (*,*) "alpha=", alpha
!stop
!---------------------------------------------------

            do i=1,N
              T_new(i)=T_new(i) + alpha*h(i)
            enddo

! write (*,*) t_new
!stop
!-------------------------  g_new  -----------------

  g_new(1) = T_new(1) -T_old(1)
       do i = 2, N-1
        g_new(i) = -aa* T_new(i-1) +(2*aa + 1)*T_new(i) - aa* T_new(i+1) - T_old(i)
      enddo
   g_new(N) = T_new(N) - T_old(N) 
 !write (*,*) G_new
!stop
!--------------------------------------------
        beta= dot_product(g_new,g_new)/dot_product(g,g)
          ! write(*,*) "beta = ", beta
!----------------------------------------------
          do i=1,N
            h(i)= -g_new(i) + beta*h(i)
          enddo
!------------------------------------------------
          gnrm=sqrt(sum(g_new**2))
           write (*,*) it, gnrm
           if (gnrm < 1.d-7) exit  
           g=g_new


             
    enddo   !!!!!!!!it

 !stop
    !do j = 1, N
       !write(66+time,*)  T_new(j)
    !enddo

              T_old= T_new

 enddo !!!!!!!!!!time

end program heat_CG



