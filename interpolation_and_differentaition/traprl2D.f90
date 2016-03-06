program test_traprl
! 2 Dimensional Composite Tramzoidal Rule
implicit none
real(8) :: a,b,c,d,s
integer :: M,N

 a = 0.d0
 b = 1.d0
 M = 100
 c = 0.d0 
 d = 1.d0
 N = 100 
call traprl2D(f,a,b,M,c,d,N,s)
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
subroutine traprl2D(f,a,b,M,c,d,N,s)
! To approximate the integral f(x,y) 
! x over [a,b] y over [c,d]
! using composite Tramzoidal Rule
! x(i) = x(0) +i*h i=0,...M, x(0) = a , x(M)= b , h=(b-a)/M
! y(i) = y(0) +j*k j=0,...N, y(0) = c , y(N)= d , k=(d-c)/N

implicit none
real(8) :: h,k,a,b,c,d,x,y,f
real(8) :: S,s1,s2,s3,s4,s5
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

do i = 1, (M-1)
    x = a + i * h
    s1 = s1 + f(x,c)
    s2 = s2 + f(x,d)
enddo

do j = 1,(N-1)
    y = c + j* k
    s3 = s3 + f(a,y)
    s4 = s4 + f(b,y)
enddo
do j = 1,(M-1)
    do i = 1,(N-1)
        x = a + i * h 
        y = c + j * k
        s5 = s5 + f(x,y)
    enddo
enddo

S = 25.d-2 *h*k*(f(a,c)+f(b,c)+f(a,d)+f(b,d) + 2.d0 &
    * ( s1 + s2 + s3 + s4) + 4.d0 * s5 ) 

end subroutine
!--------------------------------------
end program