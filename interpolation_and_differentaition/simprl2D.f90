program test_simprl
! 2 Dimensional Composite Simpson's Rule !!!Wrong answer
implicit none
real(8) :: a,b,c,d,s
integer :: M,N

 a = 0.d0
 b = 1.d0
 M = 100
 c = 0.d0 
 d = 1.d0
 N = 100 
call simprl2D(f,a,b,M,c,d,N,s)
write(*,'("Solution :" ,f22.14)') S

 contains
!-------------------------------------- 
function f(x,y)
real(8),intent(in) :: x,y
real(8) :: f,x2,y2,y4
x2 = x*x
y2 = y*y
y4 = y2*y2
f = 8.d0 / exp(x2 + y4)
end function
!--------------------------------------
subroutine simprl2D(f,a,b,M,c,d,N,s)
! To approximate the integral f(x,y) 
! x over [a,b] y over [c,d]
! using composite Simpson Rule
! x(i) = x(0) +i*h i=0,...M, x(0) = a , x(M)= b , h=(b-a)/M
! y(i) = y(0) +j*k j=0,...N, y(0) = c , y(N)= d , k=(d-c)/N

implicit none
real(8) :: h,k,a,b,c,d,x,y,f
real(8) :: S,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12
integer :: M,N,i,j

h = (b-a)/real(M)
k = (d-c)/real(N)
write(*,110) h,k
110 format ("d_x :",f8.5,3x," d_y : ",f8.5)
s = 0.d0
s1 = 0.d0
s2 = 0.d0
s3 = 0.d0
s4 = 0.d0
s5 = 0.d0
s6 = 0.d0
s7 = 0.d0
s8 = 0.d0
s9 = 0.d0
s10 = 0.d0
s11 = 0.d0
s12 = 0.d0
do j = 1,N
   y = c + (2*j-1)*k
   s1 = s1 + f(a,y)
   s3 = s3 + f(b,y)
enddo
do j = 1,N-1
   y = c + 2*j*k
   s2 = s2 + f(a,y)
   s4 = s4 + f(b,y)
enddo
do i = 1,M
   x = a + (2*i-1)*h
   s5 = s5 + f(x,c)
   s7 = s7 + f(x,d)
enddo
do i = 1,M-1
   x = x + 2*i*h
   s6 = s6 + f(x,c)
   s8 = s8 + f(x,d)
enddo
do j = 1,N
   do i = 1,M 
      x=a+(2*i-1)*h
      y=c+(2*j-1)*k
      s9=s9+f(x,y)
   enddo
enddo
do j = 1,N-1
   do i = 1,M 
      x = a + ( 2 * i - 1) * h
      y = c + 2 * j * k
      s10 = s10 + f(x,y)
   enddo
enddo

do j = 1,N
   do i = 1,M-1
      x=a+(2*i)*h
      y=c+(2*j-1)*k
      s11=s11+f(x,y)
   enddo
enddo

do j = 1,N-1
   do i = 1,M-1
      x=a+(2*i)*h
      y=c+(2*j)*k
      s12=s12+f(x,y)
   enddo
enddo


S = 1.d0/9.d0*h*k*(f(a,c)+f(b,c)+f(a,d)+f(b,d) + 4.d0* &
     s1 + 2.d0*s2 + 4.d0 *s3 + 2.d0*s4 + 4.d0 * s5 +&
     2.d0*s6+4.d0*s7+2.d0*s8+16.d0*s9+ 8.d0*s10+ 8.d0*s11+ 4.d0*s12) 

end subroutine
!--------------------------------------
end program