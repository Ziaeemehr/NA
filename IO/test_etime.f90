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
