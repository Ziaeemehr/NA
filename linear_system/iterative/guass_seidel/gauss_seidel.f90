program main
implicit none
!    Input, real(8) A(N,N), the matrix.
!    Input, real(8) B(N), the right hand side.
!    Input, real(8) X(N), the current solution estimate.
!    with a small modification of the jacobi method x(i) is used 
!    immediately after it has been evaluated. not suitable for 
!    parallel computations
integer :: i,j,p,it, it_max
integer ,parameter :: n = 4
real(8) :: r_max, d_max,tol, s, norm
real(8) :: A(n,n), b(n), x(n), x_new(n), r(n)
tol = 1.D-11
! Matrix A
A(1,1) = 7.;  A(1,2) = 3.;  A(1,3)= -1.; A(1,4)=2.
A(2,1) = 3.;  A(2,2) = 8.;  A(2,3)= 1.;  A(2,4)=-4.
A(3,1) = -1.; A(3,2) =1.;   A(3,3)= 4.;  A(3,4)=-1.
A(4,1) = 2.;  A(4,2) =-4.;  A(4,3)= -1.; A(4,4)=6.

b = (/-1., 0., -3., 1./)
x = (/0., 0., 0.,0./)
write(*,*) ' Solving Ax=b problem with Gauss_Seidel method '
write(*,*) "A = "
do i=1,n
    write (*,201) (A(i,j),j=1,n)
end do
201 format (4f12.6)

it_max = 200

do it=1,it_max
    do i = 1,n
        s =0.
        do j=1,n
            if (i==j) then
                cycle
            else
                s = s + A(i,j) *x(j)
            endif
        enddo
        x(i) = (b(i)-s)/A(i,i)
    enddo
    norm =0.
    r = matmul ( a(1:n,1:n), x(1:n) )
    r(1:n) = b(1:n) - r(1:n)
    r_max = maxval ( abs ( r ) )
    write ( *, '(2x,i6,2x,g14.6)' ) it, r_max
    if (r_max < tol) then
    write ( *, '(2x,A6,2x,A7)') "it" , "r_max"
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Convergence criterion satisifed on step ', it
        exit
    endif
enddo
write (*, '(A4)') "X = "
write (*,201) (x(j),j=1,n)
end program
