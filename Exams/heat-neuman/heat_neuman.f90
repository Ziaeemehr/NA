program heat_rod
!
! To solve the heat eauation for a rod with 
! Tt(x,t) =sigma* Txx(x,t)
! T(x,t) = 0 and T(L/4,t)= 100, 
! N.B.C : Tx(0,t) = Tx(L,t) = 0
! x  over [0,1] and t over [0,0.1]
! Variables :
! dl, dd, du subdiagonal, maindiagonal and superdiagonal elements respectively
! To plot the result:
! echo 'splot "data" nonuniform matrix with lines' | gnuplot --persist

implicit none

integer :: n,m         ! n=11, m=31
integer :: i,j,index_x
real(8) :: alph,l,time,d_x,d_t,sigma
real(8), allocatable :: dl(:),dd(:),du(:)
real(8), allocatable :: T_new(:),T_old(:)


open(unit=7, file="T.txt")

l = 1.d0
!time = 1.d-1
!d_x = l/real(n-1)
d_x = 5.d-2
n = l/d_x + 1
m = 200
!d_t  = time/real(m-1)
d_t = 125.d-5
!print *, n,m
allocate (T_new(n),T_old(n),dl(n),du(n),dd(n))

sigma = 1.d0 ! coefficient of termal conductivity
alph = (sigma * d_t)/(d_x*d_x)
!print *, alph
if (alph > 0.5) then
print*, "Solution may be unstable"
!stop
endif

du = -alph
dl = -alph
index_x = l/4 *n
!print *, index_x
!du(1) = 0.d0         !Dirichlet.B.C
du(1) = -2.d0 * alph  !Neumann B.C
!dl(n) = 0.d0         !Dirichlet B.C
dl(n) = -2.d0 * alph  !Neumann B.C

dd = 1.d0 + 2.d0 * alph
dd(index_x) = 1.d0
dl(index_x) = 0.d0
du(index_x) = 0.d0
!dd(1) = 1.d0
!dd(n) = 1.d0
!dd(1) = 1.d0 + 2.d0 * alph
!dd(n) = 1.d0 + 2.d0 * alph

T_old = 0.d0

!print *, index_x
T_old(index_x) = 100.
!print *,T_old

! print the x positions at the first row
write(7,120) n,((i-1)*d_x,i=1,n)
120 format (I3,999f12.2)

write(7,110) 0.,(T_old(i),i=1,n)


do j = 1,m
   call tridag(dl,dd,du,T_old,T_new,n)
   write(7,110) j*d_t ,(T_new(i),i=1,n)
   T_old = T_new
enddo

! Printing the results
!write(*,120) ((i-1)*h,i=2,nx-1)
!120 format (7x,' t',7x,'x = ',f5.2,9f11.2)
!do j = 1,m
!    write(7,110) (T_plot(i,j),i=1,n)
!enddo
110 format (f6.3,999f12.6)

close(7)

 contains
 
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
do  j=2,n
!Decomposition and forward substitution.
gam(j)=c(j-1)/bet
bet=b(j)-a(j)*gam(j)
if(bet.eq.0.)pause 'tridag failed'
!Algorithm fails; see below.
u(j)=(r(j)-a(j)*u(j-1))/bet
enddo 
do  j=n-1,1,-1
!Backsubstitution.
    u(j)=u(j)-gam(j+1)*u(j+1)
enddo 
return
END SUBROUTINE
end program