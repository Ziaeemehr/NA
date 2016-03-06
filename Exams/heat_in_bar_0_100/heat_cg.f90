pROGRAM heat_rod
!
! To solve the heat eauation for a rod with 
! Tt(x,t) =sigma* Txx(x,t)
! T(x,t) = 0 and T(0,t)= 100, 
! x  over [0,1] and t over [0,.1]
! Variables :
! dl, dd, du subdiagonal, maindiagonal and superdiagonal elements respectively
! To plot the result:
! echo 'splot "T.txt" nonuniform matrix with lines' | gnuplot --persist

IMPLICIT NONE

integer :: N
integer :: i,j,it,t_iter
real(8) :: alph,l,beta,denom,a,s,tol
real(8) :: gnorm,gg
real(8) :: dx,dt,sigma
real(8),allocatable :: T_new(:),T_old(:),g_new(:),g(:),h(:)

N = 30

allocate (T_new(N),T_old(N),g_new(N),g(N),h(N))
l = 1.d0
dx = l/dble(N-1)
a = 5.d-1
sigma = 0.2
dt = a * dx*dx/sigma
s = 1.d0 + 2.d0 * a
tol = 1.d-7
!print *, dt

! Initial Conditions
T_old = 0.d0
T_old(1) = 100.d0
T_new = 0.d0
T_new(1) = 100.d0
T_new(N) = 0.d0
open(unit=7, file="cg.txt")
! print the x positions at the first row
write(7,110) N,((i-1)*dx,i=1,N)
110 format (I14,*(f14.7))

DO t_iter = 1,100                    !time iteration
    g(1) = T_new(1) - T_old(1)
    DO i = 2,N-1
        g(i) = - a * T_new(i+1) + s * T_new(i)- a * T_new(i-1) - T_old(i)
    ENDDO
    g(N) = T_new(N) - T_old(N)
    ! CG LOOP
    h = - g    
    DO it = 1, 10000                ! cg iteration
        denom = -a * (h(1) * h(2) - h(N-1) * h(N)) + h(1)*h(1) + h(N)*h(N)
        DO i = 2,N-1
            denom = denom + -2.d0 * a * (h(i) * h(i+1)) + &
            s * h(i) * h(i)
        ENDDO
        gg = DOT_PRODUCT(g,g)
        alph = gg / denom
        
        T_new = T_new + alph * h
        
        ! update g
        g_new(1) = g(1) + alph*h(1)
        do i=2,N-1
           g_new(i)=g(i)+alph*(-a*h(i-1)+s*h(i)-a*h(i+1))
        enddo
        g_new(N)=g(N)+alph*h(N)
        
        beta = DOT_PRODUCT(g_new,g_new)/gg        
        h = -g_new + beta * h
        gnorm = SQRT(DOT_PRODUCT(g_new,g_new))
        print *, it, gnorm
        IF (gnorm < tol) EXIT        
        g = g_new        
    ENDDO
    
    T_old = T_new
    write(7,120) (t_iter-1)*dt , (T_new(i),i=1,N)
ENDDO

! print *, "          it          g_norm"
120 format (f14.5,*(f14.7))
close(7)

END PROGRAM
