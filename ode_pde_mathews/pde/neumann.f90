program neumann
! To approximate the laplace equation in 2D
implicit none
integer  :: n,m
real(8), allocatable :: U(:,:)
real(8)  :: a,b,h,ave,x
real(8)  :: relax,err,w,pi,tol
integer  :: i,j,cnt,max_it
pi = 4.d0 * atan(1.d0)
a = 4.d0
b = 4.d0
h = 5.d-1
n = a/h + 1
m = b/h + 1
allocate(U(n,m))
tol = 1.d-3
max_it = 100
ave = (a *(f1(0.d0) + f2(0.d0)) + b * (f3(0.d0) + f4(0.d0)))/(2.*a + 2.*b)
U(1:n,1:m) = ave
! B.C
do i = 1,m
    x = (i-1) * h
    U(1,i) = f3(x)
    U(n,i) = f4(x)
enddo

do i = 1,n
    x = (i-1) * x
    !U(i,1) = f1(x)
    U(i,m) = f2(x)
enddo

!U(1,1) = (U(1,2)   + U(2,1))  * 5.d-1
U(1,m) = (U(1,m-1) + U(2,m))  * 5.d-1
!U(n,1) = (U(n-1,1) + U(n,2))  * 5.d-1
U(n,m) = (U(n-1,m) + U(n,m-1))* 5.d-1

! SOR parameter
w = 4.d0/(2.d0 + sqrt(4.d0-(cos(pi/(n-1))+cos(pi/(m-1)))**2))
err = 1.d0
 cnt = 0
do while ((err> tol).and.(cnt<=max_it))
   err = 0.d0
   do j = 1,m-1
       do i = 2,n-1
           if (j==1) then
               relax = w * (2.d0 * U(i,j+1) + U(i+1,j) + U(i-1,j)- 4. *U(i,j)) * 0.25
               U(i,j) = U(i,j) + relax
           else
               relax = w * (U(i,j+1) + U(i,j-1) + U(i+1,j) + U(i-1,j)- 4. * U(i,j)) * 0.25
               U(i,j) = U(i,j) + relax
           endif
           if (err <= abs(relax)) then
               err = abs(relax)
           endif
       enddo
   enddo
 cnt = cnt +1
!print *, err
end do
print *, cnt    ! print the number of iterations
!write(*,120) ((i-1)*h,i=2,nx-1)
!120 format (7x,' t',7x,'x = ',f5.2,9f11.2)
do j = m,1,-1
    write(*,110) (U(i,j),i=1,n)
enddo
110 format (999f12.4)
deallocate(U)
! U = flipud(U');

 contains
 function f1(x)
     real(8),intent(in) :: x
     real(8) :: f1
     f1 = 2.d1 
 end function f1
 function f2(x)
     real(8),intent(in) :: x
     real(8) :: f2
     f2 = 18.d1
 end function f2
 
 function f3(x)
     real(8),intent(in) :: x
     real(8) :: f3
     f3 = 8.d1
 end function f3
 
 function f4(x)
     real(8),intent(in) :: x
     real(8) :: f4
     f4 = 0.d0
 end function f4
 
end program


       
           
    

