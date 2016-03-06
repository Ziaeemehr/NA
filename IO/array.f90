program main
! To test the accessibility of the columns and rows in fortran
!       Abolfazl Ziaeemehr 
!
!       Feb, 2016, IASBS
!
implicit none
integer :: i,j,k
integer, parameter :: n = 600
integer :: a(n,n,n)
integer :: clck_counts_beg, clck_counts_end, clck_rate
real(8) :: elapsed_time
call system_clock ( clck_counts_beg, clck_rate )

do i = 1,n
    do j = 1,n
        do k = 1,n
            a(i,j,k) = i + j + k 
        enddo
    enddo
enddo

call system_clock ( clck_counts_end, clck_rate )
elapsed_time = (clck_counts_end - clck_counts_beg) / real (clck_rate)
write(*,100) "i  j  k" , elapsed_time


call system_clock ( clck_counts_beg, clck_rate )

do i = 1,n
    do k = 1,n
        do j = 1,n
            a(i,j,k) = i + j + k 
        enddo
    enddo
enddo


call system_clock ( clck_counts_end, clck_rate )
elapsed_time = (clck_counts_end - clck_counts_beg) / real (clck_rate)
write(*,100) "i  k  j" , elapsed_time


call system_clock ( clck_counts_beg, clck_rate )

do j = 1,n
    do i = 1,n
        do k = 1,n
            a(i,j,k) = i + j + k 
        enddo
    enddo
enddo


call system_clock ( clck_counts_end, clck_rate )
elapsed_time = (clck_counts_end - clck_counts_beg) / real (clck_rate)
write(*,100) "j  i  k" , elapsed_time


call system_clock ( clck_counts_beg, clck_rate )

do j = 1,n
    do k = 1,n
        do i = 1,n
            a(i,j,k) = i + j + k 
        enddo
    enddo
enddo


call system_clock ( clck_counts_end, clck_rate )
elapsed_time = (clck_counts_end - clck_counts_beg) / real (clck_rate)
write(*,100) "j  k  i" , elapsed_time


call system_clock ( clck_counts_beg, clck_rate )

do k = 1,n
    do i = 1,n
        do j = 1,n
            a(i,j,k) = i + j + k 
        enddo
    enddo
enddo


call system_clock ( clck_counts_end, clck_rate )
elapsed_time = (clck_counts_end - clck_counts_beg) / real (clck_rate)
write(*,100) "k  i  j" , elapsed_time


call system_clock ( clck_counts_beg, clck_rate )

do k = 1,n
    do j = 1,n
        do i = 1,n
            a(i,j,k) = i + j + k 
        enddo
    enddo
enddo


call system_clock ( clck_counts_end, clck_rate )
elapsed_time = (clck_counts_end - clck_counts_beg) / real (clck_rate)
write(*,100) "k  j  i" , elapsed_time

100 format (A20,f12.6)

end program

! reaulst for n = 600
!             i  j  k   10.734000
!             i  k  j    9.327000
!             j  i  k    3.299000
!             j  k  i    1.774000
!             k  i  j    3.354000
!             k  j  i    1.767000




