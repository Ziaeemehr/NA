program test_linsht
! 
! Linear Shooting method
! To approximate the solution of boundary value problem 
!          x" = p(t)x'(t)+q(t)x(t)+r(t)
!          x(b) = beta in [a,b] 
! and using Runge-Kutta method of order N = 4
! 
! to plot the result
! plot "result" u 2:3 title "U", "result" u 2:4 title 'W','result' u 2:5 title "X"

implicit none
    integer,parameter :: n = 2
    integer,parameter :: nstep = 20
    real(kind=8) :: U(nstep)             ! U contains the solutions at each step
    real(kind=8) :: V(nstep),out(nstep)  ! V just like U for another equation 
    real(kind=8) :: a,b,h,alph,beta,w
    real(kind=8) :: x(0:n)
    integer :: i
    a = 0.d0
    b = 4.0
    alph = 125.d-2
    beta  = -95.d-2
    h = (b-a)/nstep
    ! x = t0 x0 x'0
    open (unit = 7, file = "result_linsht")
    x = (/0.d0, alph, 0.d0/)

    call rks4(n,h,x,nstep,xprsysU,U)
    a = 0.d0
    b = 4.0
    x = (/0.d0, 0.d0, 1.d0 /)
    h = (b-a)/nstep

    call rks4(n,h,x,nstep,xprsysV,V)

    ! calculate the solution to the boundary value broblem
    w = (beta - U(nstep))/V(nstep)
    !write (*,*) w,U(nstep),V(nstep)
    !w = 0.485884
    out = U + w * V
    do i = 1,nstep
        write(7,100) i,i*h,U(i), w*V(i), out(i)
    enddo
    100 format (1x,I5,f10.2,f15.6,f15.6,f15.6)
    close(7)

    contains

subroutine xprsysU(n,x,f)
! x prime of system of equations

    real(8), dimension(0:n) :: x,f
    integer :: n
    ! time is introduced as a new differential equation
    f(0) = 1.d0
    f(1) = x(2)
    f(2) = 2.*x(0)/(1+x(0)*x(0))*x(2)- 2./(1.+x(0)*x(0))*x(1) + 1.0
end subroutine

subroutine xprsysV(n,x,f)
! x prime of system of equations

    real(kind=8), dimension(0:n) :: x,f
    integer :: n
    ! time is introduced as a new differential equation
    f(0) = 1.d0
    f(1) = x(2)
    f(2) = 2.*x(0)/(1+x(0)*x(0))*x(2)- 2./(1.+x(0)*x(0))*x(1) 
end subroutine

subroutine rks4(n,h,x,nstep,sub,xout)
implicit none
    real(kind=8) :: x(0:n)
    real(kind=8) :: y(0:n), f(0:n,4)  
    real(kind=8) :: xout(nstep) 
    integer :: i,k,n,nstep
    real(kind=8) :: h      
    interface
        subroutine sub(n,x,f)
            real(8), dimension(0:n) :: x,f
            integer :: n
        end subroutine sub 
    end interface      
    do k = 1,nstep
        call sub(n,x,f(0,1))
        y(:) = x(:) + 0.5*h*f(:,1)
        call sub(n,y,f(0:n,2))
        y(:) = x(:) + 0.5*h*f(:,2)
        call sub(n,y,f(0:n,3))    
        y(:) = x(:) + h*f(:,3)
        call sub(n,y,f(0:n,4))    
        x(:) = x(:) + (h/6.0)* (f(:,1) + 2.0*(f(:,2) + f(:,3)) + f(:,4)) 
        xout(k) = x(1)          
    end do       
end subroutine rks4

end program