program test
! to test the random module in random.f90
use random
implicit none
real:: a
integer :: i,seed
seed = 123456
open (unit = 21, file = "result")
do i = 1,10000
  !a =  random_normal()
  a= random_gamma(1.0, .false.)
  write(21,110)  a
enddo
110 format (f10.7)
end program
