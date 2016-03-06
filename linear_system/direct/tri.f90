program main      
! to solve a tridiagonal matrix
!
!    Abolfazl Ziaeemehr
!
!    Date : 06.02.2016
!
implicit none
integer,parameter ::n = 31
real, dimension (n):: a,d,c,b,x
integer :: i
      
do i=1,n 
   d(i) = 2.0    
   a(i) = 1
   c(i) = 1 
end do
b = (/1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,15.,14.,13.,12.,11.,10.,9.,8.,7.,6.,5.,4.,3.,2.,1./)
    
call tri(n,a,d,c,b,x) 
print*, "subroutine tri"
do i = 1,n
    write(*,100) x(i)  
enddo

100 format (f12.5)
  
end program main

subroutine tri(n,a,d,c,b,x) 
!
!     a is subdiagonal of the coefficient matrix
!     d is main diagonal of the coefficient matrix
!     c is superdiagonal of the coefficient matrix
!     b is the constant vector of the linear system
!     x is the solution vector

implicit none
integer, intent(in)::n                                       
real, dimension(n) ::a,d,c,b                      
real, dimension(n) :: x                          
integer :: i                                                  
real :: mult
do i = 2,n
    mult = a(i-1)/d(i-1)
    d(i) = d(i) - mult * c(i-1)
    b(i) = b(i) - mult * b(i-1)
enddo
x(n) = b(n)/ d(n)
do i = n-1,1,-1
    x(i) =  (b(i)-c(i) * x(i+1))/d(i)
enddo
end subroutine    