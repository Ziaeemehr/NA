program test_rks4
! for table 9.14 mathews 3rd Ed page 491
implicit none
integer,parameter :: n = 2
integer,parameter :: nstep = 50
real(8) :: a,b,h,t,y
real(8) :: x(0:n),xout(nstep,0:1)
integer :: i
a = 0.d0
b = 5.d0
h = (b-a)/nstep
x = (/0.0, 3.0, -5.0/)

call rks4(n,h,x,nstep,xout)

do i = 1,nstep
    t = a + h * i
    y = 3. /exp(2.*t) * cos(t) + 1./exp(2*t) * sin(t)
    write(*,300) i,xout(i,0), xout(i,1), y
    300 format (1x,I5,f10.2,f15.8,f15.8)
enddo
end program

subroutine xprsys(n,x,f)
! x prime of system of equations
! for table 9.13 mathews 3rd Ed page 489
    real(8), dimension(0:n) :: x,f
    integer :: n
    ! time is introduced as a new differential equation
    f(0) = 1.d0
    f(1) = x(2)
    f(2) = -5.* x(1) - 4. * x(2)
end subroutine

subroutine rks4(n,h,x,nstep,xout)
! 
! n            the number of differential equations
! h            time step
! xout         x array at time tk
implicit none
      real(8) :: x(0:n)
      real(8) :: y(0:n), f(0:n,4)
      real(8) :: xout(nstep,0:1) 
      integer :: i,k,n,nstep
      real(8) :: h
      write(*,200) "k","t","x","y"
      write(*,100) 0,x(0),x(1)
      f = 0.
out:  do k = 1,nstep
          call xprsys(n,x,f(0,1))
in1:      do i = 0,n
              y(i) = x(i) + 0.5*h*f(i,1)
              !print *, f(i,1)
          end do in1          
          call xprsys(n,y,f(0:n,2))
in2:      do i = 0,n
              y(i) = x(i) + 0.5*h*f(i,2)      
          end do in2
          call xprsys(n,y,f(0:n,3))    
in3:      do i = 0,n
              y(i) = x(i) + h*f(i,3)       
          end do in3
          call xprsys(n,y,f(0:n,4))    
in4:      do i = 0,n
              x(i) = x(i) + (h/6.0)* (f(i,1) + 2.0*(f(i,2) + f(i,3)) + f(i,4)) 
          end do in4
          xout(k,0:1) = x(0:1)
          !write(*,100) k,x(0),x(1)
      end do out
      100 format (1x,I5,f10.2,f15.8,f15.8)
      200 format (1x,A4,2A10,A16)
end subroutine rks4 
