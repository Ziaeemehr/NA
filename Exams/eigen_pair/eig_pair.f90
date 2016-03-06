program eigpair
!
!    Abolfazl Ziaeemehr
!
!    Date : Feb.2016
!
!    To solve Ax=b from calculation of eigen-pair of the matrix 
!    using MKL subroutines ddot and dsyev
!

implicit none
!include 'mkl.fi'
!include 'blas.f90'
integer,parameter :: N = 5
real(kind=8) :: b(N)
real(kind=8) :: ddot
real(kind=8) :: a(5,5) ! Input: The coefficient matrix / Output: The eigenvectors
real(kind=8) :: work(3*N-1)
real(kind=8) :: w(N) ! contains the eigenvalues of the matrix A in ascending order
real(kind=8) :: x(N), y(N)
character(len=1) :: jobz
character(len=1) :: uplo
integer :: lda
integer :: lwork
integer :: incx, incy,i,j,info
real(kind=8) :: res

jobz = "V"
uplo = "U"
lda = max(1,N)
lwork = max (1, 3*N-1)
a(1,1:5) = (/5.d0,-2.d0,0.d0,-1.d0,0.d0/)
a(2,1:5) = (/-2.d0,5.d0,0.d0,0.d0,-1.d0/)
a(3,1:5) = (/0.d0,0.d0,5.d0,0.d0,0.d0  /)
a(4,1:5) = (/-1.d0,0.d0,0.d0,5.d0,-2.d0/)
a(5,1:5) = (/0.d0,-1.d0,0.d0,-2.d0,5.d0/)

b(1:5) = (/1.d0,2.d0,3.d0,4.d0,5.d0/)

! print the input data
print *, "Matrix A is "
do i =1,5
    write(*,100) (a(i,j), j=1,5)
end do

print *, " b is "
write (*,'(f10.5)') (b(i), i = 1,N)
100 format (5f10.2)

call dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)

print *, "The eigenvalues in ascending order "
write(*,110) (w(i),i=1,N)
110 format (5f12.6)

print *, "The eigenvectors are : "
do i = 1,N
    write(*,100) (a(i,j), j=1,5)
enddo

incx = 1
incy = 1

do i = 1,N
    res = ddot (N,a(:,i),incx,b(:),incy)
    x = x + 1.d0/w(i) * res * a(:,i)
enddo
print *, ""
print *, "The x vector is : "
write (*,'(f10.5)') (x(i),i=1,N)

end program
