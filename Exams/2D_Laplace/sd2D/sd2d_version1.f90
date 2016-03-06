PROGRAM SD
!
! purpose:
! Solving 1D poisson equation from Steepest Descent method
! Ploting the results
! echo 'splot "v.txt" nonuniform matrix w l lt 8 title "SD"' | gnuplot --persist
IMPLICIT NONE
REAL(8), allocatable :: v(:,:), g(:,:),v_ex(:,:)
INTEGER :: n,m,i,j,iter,ii,jj
REAL(8) :: pi,alph,h,k,x,y,norm_g,tol,a,b
INTEGER :: ieror
pi = 4.d0 * atan(1.d0)
! number of intervals

a = 100.d0
b = 1.d0 

m = 1000
tol = 1.d-2
h = 1.d2/real(m)
n=m/100
!k = 1.d0/real(n)
!print *, h
!stop
allocate (v(m+1,n+1),g(m+1,n+1),V_ex(m+1,n+1))
! Boundary Conditions:
v = 0.d0
do i = 1,m+1
    x = (i-1)*h
    v (i,n+1) = sin(pi*x/a)
enddo
!print *,maxval(v(:,n+1))
!stop
alph = 25.d-2     
! h    = 1.d0/n

DO i = 1,m+1
    DO j = 1,n+1
        x = (i-1) * h
        y = (j-1) * h
        v_ex(i,j) = sin(pi*x/a)*sinh(pi*y/a)/sinh(pi*b/a)
    ENDDO
ENDDO
! print some variables to check
print *, "y=", y,"h=",h,"n=", n

DO iter  = 1, 1000000    
    g(1,1) = 4.d0 * v(2,2) -v(3,2) - v(1,2) - v(2,3) - v(2,1)
    DO i = 3,m
        DO j = 3,n    
            g(i-1,j-1) = 4.d0 * v(i,j) -v(i+1,j)-v(i-1,j)-v(i,j+1)-v(i,j-1)
        ENDDO
    ENDDO    
    g(m,n) = 4.d0 * v(m,n) -v(m-1,n)-v(m+1,n)-v(m,n+1)-v(m,n-1)
    norm_g = 0.d0
    DO ii = 1,m-1
        DO jj = 1,n-1
           norm_g = norm_g + g(ii,jj) * g(ii,jj)
        ENDDO
    ENDDO
    norm_g = sqrt(norm_g)
    IF (mod(iter,100)==0)   print *, iter,norm_g
    !print *, iter, norm_g
    IF(norm_g < tol) exit
    v(2:m,2:n) = v(2:m,2:n) - alph * g
ENDDO
!print *,'	It  	Norm of the gradient'
write(*,*) iter, norm_g
open (8, FILE = 'v.txt', iostat = ieror)
write(8,100) N+1,(i*h,i = 0,n)
do i = 1, m+1
    write(8,110) (i-1)*h, (v(i,j),j = 1,n+1)
enddo
100 format (I5,*(f12.7))
110 format (f10.4,*(f12.7))
! To check the difference of numerical and analytical solutions
print *, "V - V_exact"
do j = 1,n+1
print *, maxval(v(:,j)-v_ex(:,j))
enddo
close(8)
END PROGRAM SD
