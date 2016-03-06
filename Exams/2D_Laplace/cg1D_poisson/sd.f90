PROGRAM SD
!
! purpose:
! Solving 1D poisson equation from Steepest Descent method
! Ploting the results
! echo 'plot "v1.txt" u 1:2 w l title "SD" ,"v1.txt" u 1:3 w l title "EXACT"' | gp --persist
IMPLICIT NONE
REAL(8), allocatable :: v(:), g(:), rho(:),V_ex(:)
INTEGER :: n, ip, ii, it
REAL(8) :: pi , alph, h, x, norm_g,tol
INTEGER :: ieror
! number of intervals
n = 200
tol = 1.d-9
allocate (v(n-1),g(n-1),rho(n-1),V_ex(n-1))

!      print *
!      print *,' Steepest Descent method'
!      print *

! Boundary Conditions:
!      v(0) = 0
!      v(n) = 1

! Steepest descent length scale
alph = 25.d-2     
h    = 1.d0/n
pi   = 4.d0*atan(1.d0)

DO ip = 1,n-1
    x = ip * h
    rho(ip) = 10.d0 * sin(10.d0 * sin(pi*x)) * (cos(pi*x))**2 + sin(pi*x) * cos(10.d0 * sin(pi*x))
    V_ex(ip) = 1.d-1*(sin(10.d0 * sin(pi*x)))/pi**2+x
ENDDO

DO ii  = 1, 1000000
    g(1) = 2 * v(1) - v(2) - h**2 * rho(1)
    DO it  = 2,n-2
        g(it) = -v(it-1)+2*v(it)-v(it+1)-h**2*rho(it)
    ENDDO
    g(n-1) = -v(n-2) + 2*v(n-1) -1 - h**2 * rho(n-1)
    norm_g = 0.d0
    DO   ip = 1,n-1
         norm_g = norm_g + g(ip) * g(ip)
    ENDDO
    !norm_g = sqrt(sum(g(1:n-1) * g(1:n-1)))
    IF(norm_g < tol) exit
    v(1:n-1) = v(1:n-1) - alph * g(1:n-1)  
ENDDO
!print *,'	It  	Norm of the gradient'
!write(*,*) ii, norm_g
open (8, FILE = 'v1.txt', iostat = ieror)
do ii = 1,n-1
    write(8,110) ii*h, v(ii), V_ex(ii)
enddo
110 format (f10.5,2f14.9)
close(8)
END PROGRAM SD
