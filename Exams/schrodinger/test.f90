program test
implicit none
complex(8) :: z
real(8) :: x,y
x = 1.d0
y = 2.d0
z = cmplx(x,y,8)

print *, z
print*, ""
print *, real(z)
print *, aimag(z)
print*, ""
print *,z*conjg(z)


end program