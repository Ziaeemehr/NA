program main
implicit none
!    Input, real(8) A(N,N), the matrix.
!    Input, real(8) B(N), the right hand side.
!    Input, real(8) X(N), the current solution estimate.
integer :: i,j,it, it_max
integer ,parameter :: n = 4
real(8) :: r_max, d_max, tol, pAp, rsnew, rsold, alpha
real(8) :: A(n,n), b(n), x(n), r(n), p(n), Ap(n)
tol = 1.D-10
! Matrix A
A(1,1) = 7.;  A(1,2) = 3.;  A(1,3)= -1.; A(1,4)=2.
A(2,1) = 3.;  A(2,2) = 8.;  A(2,3)= 1.;  A(2,4)=-4.
A(3,1) = -1.; A(3,2) = 1.;  A(3,3)= 4.;  A(3,4)=-1.
A(4,1) = 2.;  A(4,2) =-4.;  A(4,3)= -1.; A(4,4)=6.

b = (/-1., 0., -3., 1./)
x = (/0., 0., 0.,0./)
! printing the header and matrix
write(*,*) 'Solving Ax=b problem with Conjugate gradient method'
write(*,*) "A = "
do i=1,n
    write (*,201) (A(i,j),j=1,n)
end do
201 format (4f12.6)

it_max = 1000

r = matmul ( a(1:n,1:n), x(1:n) )
r = b - r
p = r
rsold = 0.
do i = 1,n
    rsold = rsold + r(i)*r(i)
enddo
write ( *, '(a)' ) ' '
! The main loop of CG method
do it = 1, it_max
    Ap = matmul( a(1:n,1:n), p(1:n) )
    pAp= DOT_PRODUCT(p,Ap)
    alpha = rsold/pAp
    x = x + alpha * p
    r = r - alpha * Ap
    rsnew = DOT_PRODUCT(r,r)
    write ( *, '(2x,i6,2x,g16.6)' ) it, sqrt(rsnew)
    if (sqrt(rsnew)< tol) then
        write ( *, '(2x,A6,2x,A19)') "it" , "residual rs=b-Ax"
        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' ) '  Convergence criterion satisifed on step ', it-1        
        exit        
    endif
    p = r + (rsnew / rsold) * p
    rsold = rsnew
    
enddo
write (*, '(A4)') "X = "
write (*,201) (x(j),j=1,n)
end program
    



!    r = matmul ( a(1:n,1:n), x(1:n) )
!    r(1:n) = b(1:n) - r(1:n)
!    r_max = maxval ( abs ( r ) )
!    write ( *, '(2x,i6,2x,g14.6)' ) it, r_max
!    if (r_max < tol) then
!    write ( *, '(2x,A6,2x,A7)') "it" , "r_max"
!    write ( *, '(a)' ) ' '
!    write ( *, '(a,i6)' ) '  Convergence criterion satisifed on step ', it
!        exit
!    endif
!enddo
!write (*, '(A4)') "X = "
!write (*,201) (x(j),j=1,n)
!end program
