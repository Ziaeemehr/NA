program test_diffnew 
!
! To approximate f'(x) numerically by constructing the Nth degree Newton polynomial
! P(x)= a0+a1(x-x0)+a2(x-x0)(x-x1)+
! a3(x-x0)(x-x1)(x-x2) + ... + aN(x-x0)...(xN-x_{N-1})
! and using f'(x0) \approx P'(x0)



subroutine diffnew(x,y,n)
implicit none
A = y
do j = 2,N
   do k = n:j,-1
       A(k) = (A(k)-A(k-1))/(x(k)-x(k-j+1))
   enddo
enddo
x0 = x(1)
df = A(2)
prod = 1
n1 = size(A)-1
do k = 2,n1
    prod = prod * (x0-x(k))
    df = df + prod * A(k+1)
end



