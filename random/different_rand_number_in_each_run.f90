program rand
!
! To test subroutine for generating different random numbers
! in each run.
implicit none
real(8):: r
real(8):: a(3)
integer :: i
call init_random_seed()
do i = 1,10
  call random_number(a)
  write(*,110) "a = ", a
enddo
100 format (A5,f10.7)
110 format (A5,3f10.7)
end program

subroutine init_random_seed()
! To produce different random numbers in each run
! seed is changed with the clock 
      INTEGER :: i, n, clock
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed

      CALL RANDOM_SEED(size = n)
      ALLOCATE(seed(n))

      CALL SYSTEM_CLOCK(COUNT=clock)

      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      CALL RANDOM_SEED(PUT = seed)

      DEALLOCATE(seed)
end
