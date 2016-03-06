program test1
implicit none
real :: result = 10.0
integer :: i = 2
character(len=20)::string
integer :: index = -12, junk = 4, number = -12345
real :: a =-12.3, b = .123, c = 123.456
real :: a1 = 1.2346E6, b1 = 0.001, c1 = -77.7E10, d1 = -77.7E10

write(*,100) i,result
100 format ('the result for iter',I5,' is',F7.3)

string = '(1X,I6,F10.2)'
write(*,string) i,result

write(*,'(1X,I6,F10.2)') i, result
!write(*,200)
!write(*,210)
!write(*,220)
!200 format ('l','This heading is at the top of a new page.')
!210 format ('O', '    Control character      Action  ')
!220 format (' ', '    =================      =====   ')


write(*,230) index, index+12, junk, number
write(*,240) index, index+12, junk, number
write(*,250) index, index+12, junk, number
230 format (' ',2I5,I6 I10)
240 format (' ',2I5.0,I6 I10.8)
250 format (' ',4I7)

write(*,*)
write(*, 260) a,b,c
write(*, 270) a,b,c
260 format (' ', 2F7.3, F8.3)
270 format (' ', 3F10.2)

! E descriptor
write(*,280) a1, b1, c1, d1
280 format (' ',2E14.4,E15.6, E11.6 )







end program

