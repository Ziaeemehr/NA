program test_adapt
! To test the adaptive simpson integration
implicit none
!interface 
!function f(x)
!    real(8),intent(in) :: x
!end function f
!end interface

real(8) ::tol,a,b,exact,err,quad
real(8) :: SRmat(30,6)

!integer ::

tol = 1.d-5
a = 0.d0
b = 4.d0
exact = -1.5487883725279481333
 
 call adapt(f,a,b,quad,SRmat,tol,err)
write(*,*) quad

 contains
!---------------------------------------------------------- 
function f(x) 
    real(8),intent(in) :: x
    real(8):: f
    f = 13.d0 *(x-x*x)*exp(-15.d-1*x)
end function f
!----------------------------------------------------------
subroutine adapt(f,a,b,quad,SRmat,tol,err)
! To approximate the integral f(x) over [a,b]
! The composite Simpson rule is applied to the 4M subintervals
! [x(4k-4),x(4k)] where [a,b]=[x(0),x(4M)] 
implicit none
integer :: n,m,iterating,j,state,p,done
real(8) :: SRmat(30,6),SRvec(6),SR0vec(6),SR1vec(6),SR2vec(6)
real(8) :: a,b,tol,f,err,quad,c,tol2

! initialize values
SRmat(30,6) = 0.d0
iterating = 0
done = 1
SRvec = 0.d0

call srule(f,a,b,tol,SRvec)
SRmat(1,1:6) = SRvec
m=1
state = iterating
do while (state == iterating)
    n = m
    do j = n,1,-1
        p = j
        SR0vec = SRmat(p,:)
        err = SR0vec(5)
        tol = SR0vec(6)
        if (tol <= err) then
	    ! Bisect interval, apply Simpson's rule
	    ! recursively, and determine error
	    state = done
	    SR1vec = SR0vec
	    SR2vec = SR0vec
	    a = SR0vec(1)
	    b = SR0vec(2)
	    c = 5.d-1 * (a+b)
	    err = SR0vec(5)
	    tol = SR0vec(6)
	    tol2 = 5.d-1 * tol
	    call srule(f,a,c,tol2,SR1vec)
	    call srule(f,c,b,tol2,SR2vec)
	    err = abs(SR0vec(3)-SR1vec(3)-SR2vec(3))/10.d0
            ! Accuracy test
	    if (err<tol) then
		SRmat(p,:) = SR0vec
		SRmat(p,4) = SR1vec(3)+SR2vec(3)
		SRmat(p,5) = err
	    else
		SRmat(p+1:m+1,:) = SRmat(p:m,:)
		m = m + 1
		SRmat(p,:)  = SR1vec
		SRmat(p+1,:)= SR2vec
		state = iterating
            endif
        endif
    enddo
enddo
quad = sum(SRmat(:,4))
err  = sum(abs(SRmat(:,5)))
SRmat = SRmat(1:m,1:6)
end subroutine
!----------------------------------------------------------
subroutine srule(f,a0,b0,tol0,z)
! To approximate the integral f(x) over [a0,b0] 
! using Simpson's rule
implicit none
real(8) :: z(6),c(3)
real(8) :: h,f,a0,b0,tol0,tol1,err,s,s2

 h = (b0 - a0) * 5.d-1
 c = 0.d0
 c(1) = f(a0)
 c(2) = f(5.d-1*(a0+b0))
 c(3) = f(b0)
 s = h*(c(1)+4.d0*c(2)+c(3))/3.d0
 s2   = s
 tol1 = tol0
 err  = tol0
 z = (/a0, b0, s, s2, err, tol1/)
 
end subroutine 
!----------------------------------------------------------
end program