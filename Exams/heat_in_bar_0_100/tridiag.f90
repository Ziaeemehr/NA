SUBROUTINE tridag(a,b,c,r,u,n)
INTEGER n,NMAX
REAL(8) :: a(n),b(n),c(n),r(n),u(n)
PARAMETER (NMAX=500)
! Solves for a vector u(1:n) of length n the tridiagonal linear set 
! a(1:n) , b(1:n) , c(1:n) , and r(1:n) are input vectors and are not modified.
! Parameter: NMAX is the maximum expected value of n .
INTEGER j
REAL(8) :: bet,gam(NMAX)      !One vector of workspace, gam is needed.
if(b(1).eq.0.)pause 'tridag: rewrite equations'
!If this happens then you should rewrite your equations as a set of order N − 1, with u2
!trivially eliminated.
bet=b(1)
u(1)=r(1)/bet
do 11 j=2,n
!Decomposition and forward substitution.
gam(j)=c(j-1)/bet
bet=b(j)-a(j)*gam(j)
if(bet.eq.0.)pause 'tridag failed'
!Algorithm fails; see below.
u(j)=(r(j)-a(j)*u(j-1))/bet
enddo 11
do 12 j=n-1,1,-1
!Backsubstitution.
u(j)=u(j)-gam(j+1)*u(j+1)
enddo 12
return
END SUBROUTINE