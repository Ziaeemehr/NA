program test_rks4
! for table 9.13 mathews 3rd Ed page 489
implicit none
integer,parameter :: n = 2
integer,parameter :: nstep = 10
real(8) :: a,b,h
real(8) :: x(0:n)

a = 0.d0
b = 2.d-1
h = (b-a)/nstep
x = (/0.0, 6.0, 4.0/)

call rks4(n,h,x,nstep)

end program

subroutine xprsys(n,x,f)
! x prime of system of equations
! for table 9.13 mathews 3rd Ed page 489
    real(8), dimension(0:n) :: x,f
    integer :: n
    ! time is introduced as a new differential equation
    f(0) = 1.d0
    f(1) = x(1) + 2 * x(2)
    f(2) = 3 * x(1) + 2 * x(2)
end subroutine


subroutine rks4(n,h,x,nstep)
! 
! To approximate the solution of the ststem of equations
! n            the number of differential equations
! h            time step
!
implicit none
      real(8) :: x(0:n)
      real(8) :: y(0:n), f(0:n,4)  
      integer :: i,k,n,nstep
      real(8) :: h
      write(*,200) "k","t","x","y"
      write(*,100) 0,x
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
          write(*,100) k,x
      end do out
      100 format (1x,I5,f10.2,f15.8,f15.8)
      200 format (1x,A4,2A10,A16)
end subroutine rks4 
