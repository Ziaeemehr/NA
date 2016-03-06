program pngauss
! to solve Ax=b by the naive gaussian elimination
integer, parameter :: n = 4 !The order of A, the number of rows in B; n≥ 0.
integer :: i,j
integer :: nrhs             ! The number of right-hand sides, the number of columns in B; nrhs≥ 0
integer :: info
real(8) :: a(n,n)
real(8) :: dl(n-1)          ! subdiagonal elements of A.
real(8) :: d(n)             ! diagonal elements of A 
                            !! Overwritten by the n diagonal elements of U as output
real(8) :: du(n-1)          ! contains the (n - 1) superdiagonal elements of A
                            !! Overwritten by the (n-1) elements of the first superdiagonal of U.

real(8) :: b(n)             ! the matrix B whose columns are the right-hand sides for the systems of equations.
                            !! Overwritten by the solution matrix X.


nrhs = 1
! initializing the matrix
!do i = 1,n
!    d(i) = 2.*i-10.
!    b(i) = i - n
!enddo
!do i = 1,n-1
!    du(i) = i*i-20.
!    dl(i) = -i*i+20.
!enddo
d = 3.
du = 4.
dl = 1.
b  = 2. 
do i = 1,n
    do j = 1,n
       if (i==j) a(i,j) =  d(i)
       if ((i+1)==j) a(i,j) = du(i)
       if (i==(j+1)) a(i,j) = dl(j)
    enddo
enddo
print *, "The coefficient matrix is: "
do i=1,n
    write(*,110) (a(i,j),j=1,n)
end do  


print*, "the following are the constants b(i)"
write(*,110) (b(i), i=1,n)


call dgtsv ( n , nrhs , dl , d , du , b , ldb , info )

print*, "the solution x(i) "

write(*,110)  (b(i), i=1,n)

110 format (1x,4f14.6)
end program
