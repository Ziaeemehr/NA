Fortran
=======
Created جمعه 12 ژوئن 2015

run:
f95 myprog.f90 -o myprog

print *, "alpha + beta = " , alpha + beta

du -sh folder will give you the size of the directory folder and du -sh . will give you the size of the current folder you are in.
ls *.f90                    # List all files whose name ends in .f90
rm *.out                  # Remove all files ending in .out
du -sh *                   # Show the size of all the files and directories
mv *.f90 programs/  # Move all files whose name ends in .f90 into programs/
ls p????.f90              # list all files whose name is p+any 4 characters+.f90
ls ?.out                    # list all files whose name is any one character followed by .out
rm -r * this would simply erase all files and folders in the current folders
 data type:
double precision :: e , f
e.g. 0d0 , -13.67d0 , 7.3d6 , 4.1d-3

Fixing Mixed Variable Types
int(a)   return the truncated value of a in integer form 13.999 → 13
nint(a)  it will return the rounded integer value so 13.999 →  14.
real(a) 13 → 13.00
dble(a) : dble(13.45) → 13.45000000000

Parameter Command
integer, parameter :: size = 20
real, parameter :: pi = 3.141
real, dimension(size) :: stars
Matrices
real, dimension(5,3) :: table

integer, parameter :: rows = 5
integer, parameter :: cols = 3
real, dimension(rows,cols) :: table

read *, table    ! we would input our values in column order
table(1,1) , table(2,1) , table(3,1) , table(4,1) , table(5,1) , table(1,2) , table(2,2) etc . . .

Implied Do-Loop

do i=1,rows
read *, (table(i,j) , j=1,cols) ! here is the implied do-loop bit
end do

Allocatable Arrays
real, allocatable, dimension(:) :: x       ! rank 1
real, allocatable, dimension(:,:) :: y     ! rank 2
print *, ’enter integer values for the size of your arrays:’
read *, a,b
allocate (x(a) , y(a,b))
deallocate(x,y)

#### External Input/Output Files
Reading in .DAT Files
open(unit=7 , file=’/home/astronomy/as3013/datfiles/numbers1.dat’ , status=’old’)
since we are in /home/astronomy/as3013  we could also use:
open(unit=7 , file=’/datfiles/numbers1.dat’ , status=’old’)
'new', if we want to create the file instead of 'old'

^{UNIT = unit number
FILE = file_name
STATUS = file_status
FORM = format_mode
ACTION = allowed_actions
POSITION = file_position
IOSTAT = ios}

Variable type
_{program test
implicit none
real :: pi1
real*8 :: pi2
real*16 :: pi3
pi1 = acos(-1.0)
pi2 = acos(-1.d0)
pi3 = acos(-1.0_16)
print *, pi1
print *, pi2
print *, pi3
end program test}

##### write
WRITE (*,100) 'The coefficients are : ', A, B, C
FORMAT (1X,A,3F10.4)

##### open
open (unit = 2, file = "scores")
open (unit = 7, file = "a:scores.txt")
open (unit = 5, file = "h:\\results\\primes")
close (2)
write (7,*)

##### read
1.23 , 4.56 , 7.89
11, 13 , "Sally"
read (7,*) x, y, z
read (7,*) m, n, first 

read (7,5) n, x, y
5 	format (i2,f5.2,f3.1)    

#### Timing in Fortran
!The following code fragment shows how to output elapsed wall clock time in seconds. This is different from consumed CPU time and may not be a good measure of code performance. DTIME is runtime and not wall clock time, also it is an extension. If you want cpu rutime the language standard provides cpu_time. Your code for converting from seconds looks OK.

integer :: clck_counts_beg, clck_counts_end, clck_rate
call system_clock ( clck_counts_beg, clck_rate )
your code
call system_clock ( clck_counts_end, clck_rate )
write (*, *)  (clck_counts_end - clck_counts_beg) / real (clck_rate)

CPU time example, also in seconds

real ::  beg_cpu_time, end_cpu_time
call cpu_time (beg_cpu_time)
your code
call cpu_time (end_cpu_time)
write (*, *) end_cpu_time - beg_cpu_time


ETIME — Execution time subroutine (or function)
syntax: 
CALL ETIME(VALUES, TIME). 
        VALUES(1): 	User time in seconds.
	VALUES(2): 	System time in seconds.
	TIME:    	Run time since start in seconds. 

program test_etime
              integer(8) :: i, j
              real, dimension(2) :: tarray
              real :: result
              call ETIME(tarray, result)
              print *, result
              print *, tarray(1)
              print *, tarray(2)
              do i=1,100000000    ! Just a delay
                  j = i * i - i
              end do
              call ETIME(tarray, result)
              print *, result
              print *, tarray(1)
              print *, tarray(2)
          end program test_etime





##### One-dimensional Arrays
A = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 /)
A = (/ (I, I = 1,10) /)

!	A scalar expression as in
	X (1:4) = (/ 1.2, 3.5, 1.1, 1.5 /)
!	An array expression as in
	X (1:4) = (/ A (I, 1:2), A (I+1, 2:3) /)
!	An implied do loop as in
	x (1:4) = (/ (sqrt (real (i)), i = 1, 4) /)

reshape ( (/ 1, 2, 3, 4, 5, 6 /), (/ 2, 3 /) )
!is the 2 x 3 array
!     1  3  5
!     2  4  6


! I/O When you don't the size of the array 
program test
    implicit none
    integer :: t(4,4)
    integer :: j
t(:, :) = 2
!open(unit=10, file='test.csv')
do j = 1, 4
   ! write(*, *) t(1:4, j)
   ! write(*, '(3(I0, ", "), I0)') t(1:4, j)
   write(*, '(*(I0 : ", "))') t(:, j)
end do
!close(10)
end program
###  Variable format
do i = 1, matrix
  write(*,'(*(i10,1x))') A
end do
###


program test_mod
            print *, mod(17,3)
            print *, mod(17.5,5.5)
            print *, mod(17.5d0,5.5)
            print *, mod(17.5,5.5d0)
            print *, ""
            print *, mod(-17,3)
            print *, mod(-17.5,5.5)
            print *, mod(-17.5d0,5.5)
            print *, mod(-17.5,5.5d0)
            print *, ""
            print *, mod(17,-3)
            print *, mod(17.5,-5.5)
            print *, mod(17.5d0,-5.5)
            print *, mod(17.5,-5.5d0)
end program test_mod



! Open input data file.  
OPEN ( UNIT=21, FILE=filename, status='OLD', ACTION='READ', &
       IOSTAT=status )
READ (21,*) nat
READ (21,*)

fileopen: IF ( status == 0 ) THEN    ! Open successful
    DO iat = 1, nat
        READ (21,*) sat, rat(1, iat), rat(2, iat), rat(3, iat)
    END DO
    CLOSE(21)   
END if fileopen



#### assumed shape array
subroutine sub90(x)
  real x(:,:)
  m=size(x,1)
  n=size(x,2)
  print *, x(m,n)
  return
end
	