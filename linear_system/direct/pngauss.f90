program pngauss
! to solve Ax=b by the naive gaussian elimination
integer, parameter :: n = 4
real(8) :: a(n,n), b(n), x(n)
integer :: i,j

a(1,1)=3.0;  a(1,2)=-13.0; a(1,3)=9.0; a(1,4)=3.0
a(2,1)=-6.0; a(2,2)=4.0;   a(2,3)=1.0; a(2,4)=-18.0
a(3,1)=6.0;  a(3,2)=-2.0;  a(3,3)=2.0; a(3,4)=4.0
a(4,1)=12.0; a(4,2)=-8.0;  a(4,3)=6.0; a(4,4)=10.0
print*, "the coefficient matrix a"
b(1)=-19.0; b(2)=-34.0; b(3)=16.0; b(4)=26.0
do i=1,n
write(*,110) (a(i,j),j=1,n)
end do
print*, "the following are the constants b(i)"
write(*,110) (b(i), i=1,n)

call ngauss(A,b,x,n)

print*, "the solution x(i) "

write(*,110)  (x(i), i=1,n)

110 format (1x,4f14.6)
end program
      

subroutine ngauss(A,b,x,n)
implicit none
integer :: i,j,k,n
real(8) :: a(n,n), b(n),x(n), xmult,sum1

x = 0.d0
do  k=1,n-1
    do  i=k+1,n      
        xmult = a(i,k)/a(k,k)       
        a(i,k) = 0.0
        do j=k+1,n    
           a(i,j) = a(i,j) - xmult * a(k,j)
        enddo
        b(i) = b(i) - xmult * b(k)
    enddo
enddo
! back substitution
x(n) = b(n)/a(n,n)
do  i = n-1,1,-1
    sum1 = b(i)
    do j = i+1,n
        sum1 = sum1 - a(i,j) * x(j)
    enddo
    x(i) = sum1 / a(i,i)
enddo

end subroutine