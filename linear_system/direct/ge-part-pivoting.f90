!!---------------------------------------------------------------------------
! Gaussian elimination with partial pivoting using straightforward formulas
! and Fortran 90/95 features
!!---------------------------------------------------------------------------
! gfortran -ffree-form -Wall -Wextra -fopenmp
!!---------------------------------------------------------------------------
! [1] J. Demmel "Numerical Linear Algebra"
! [2] N. J. Nigham "Gaussian Elimination"
!
!   ELIMINATION
!   for k = 1: n-1
!       for i = k+1: n
!           a(r(i), k) /= a(r(k), k)
!       for i = k+1: n
!           for j = k+1: n
!               a(r(i), j) -= a(r(i), k) * a(r(k), j)
!           end
!       end
!   end
!
!   BACKSUBSTITUTION  -  L U y = f  =>  L x = f  =>  x = L \ f  =>  y = U \ x
!   for i = 1: n
!       x(r(i)) = f(r(i))
!       for j = 1: i-1
!           x(r(i)) -= a(r(i), j) * x(r(j))
!       end
!   end
!
!   for i = n: 1
!       y(r(i)) = x(r(i))
!       for j = n: i+1
!           y(r(i)) -= a(r(i), j) * y(r(j))
!       end
!       y(r(i)) /= a(r(i), i)
!   end 
!!---------------------------------------------------------------------------
program ge
!!---------------------------------------------------------------------------
    use omp_lib, only: omp_get_wtime

    implicit none

    ! Matrix of coefficients; the one is filled by random_number()
    real, dimension(:, :), allocatable :: a

    ! "Analytical" solution; the one is filled by random_number()
    real, dimension(:), allocatable :: u

    ! Right-hand side (RHS); the one is calculated as f = A * u
    ! Numerical solution (NS) of the equation A y = f
    ! RHS is overwritten by NS
    real, dimension(:), allocatable :: y

    ! Order of rows
    integer, dimension(:), allocatable :: r

    ! Size of matrix
    integer, parameter :: n = 5

    ! Time stamps
    real(kind(0.d0)) :: elimination_start, elimination_finish
    real(kind(0.d0)) :: backsubstition_start, backsubstition_finish

    ! Allocate memory
    allocate(A(1: n, 1: n))
    allocate(u(1: n))
    allocate(y(1: n))
    allocate(r(1: n))

    ! Algorithm uses straightforward formulas
    call Generate_Data()
    call Straightforward_Elimination()
    call Straightforward_Backsubstition()
    call Print_Norms()
    call Print_Times()

    ! Algorithm uses Fortran 90/95 features
    call Generate_Data()
    call Fortran9x_Elimination()
    call Fortran9x_Backsubstition()
    call Print_Norms()
    call Print_Times()

    ! Free memory
    deallocate(a)
    deallocate(u)
    deallocate(y)
    deallocate(r)
!!---------------------------------------------------------------------------
contains
!!---------------------------------------------------------------------------
subroutine Print_Norms()
    write (*, *) "Norms:", maxval(abs(u)), maxval(abs(y(r) - u))
end subroutine Print_Norms
!!---------------------------------------------------------------------------
subroutine Print_Times()
    write (*, *) "Times:", &
            elimination_finish - elimination_start, &
            backsubstition_finish - backsubstition_start
end subroutine Print_Times
!!---------------------------------------------------------------------------
! This version is a simplified modification of
! http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
subroutine Init_Random_Seed()
    integer :: i, n
    integer, dimension(:), allocatable :: seed

    call random_seed(size = n)
    allocate(seed(n))

    seed = 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put = seed)

    deallocate(seed)
end subroutine Init_Random_Seed
!!---------------------------------------------------------------------------
subroutine Generate_Data()
    integer :: i

    call Init_Random_Seed()

    call random_number(A)
    call random_number(u)
    y = matmul(A, u)

    do i = 1, n
        r(i) = i
    end do
end subroutine Generate_Data
!!---------------------------------------------------------------------------
subroutine Swap_Integers(a, b)
    integer, intent(inout) :: a, b
    integer :: tmp

    tmp = a
    a = b
    b = tmp
end subroutine Swap_Integers
!!---------------------------------------------------------------------------
subroutine Straightforward_Elimination()
    integer :: i, j, k
    integer :: max_loc
    real :: max_val

    elimination_start = omp_get_wtime()

    do k = 1, n-1
        ! Find maxval(abs()) element in leading column of active submatrix
        max_loc = k
        max_val = abs(a(r(k), k))
        do i = k+1, n
            if (abs(a(r(i), k)) > max_val) then
                max_val = abs(a(r(i), k))
                max_loc = i
            end if
        end do

        ! Make that element a pivot
        call Swap_Integers(r(k), r(max_loc))

        ! Perform the rest as usual
        do i = k+1, n
            a(r(i), k) = a(r(i), k) / a(r(k), k)
        end do

        do j = k+1, n
            do i = k+1, n
                a(r(i), j) = a(r(i), j) - a(r(i), k) * a(r(k), j)
            end do
        end do
    end do

    elimination_finish = omp_get_wtime()
end subroutine Straightforward_Elimination
!!---------------------------------------------------------------------------
subroutine Fortran9x_Elimination()
    integer :: k, i
    integer :: max_loc
    real :: max_val

    elimination_start = omp_get_wtime()

    do k = 1, n-1
        ! For some reason that I don't quite understand yet, it is impossible
        ! to compile the following line:
        !max_loc = maxloc( abs( a(r(k: n), k) )  )

        max_loc = k
        max_val = abs(a(r(k), k))
        do i = k+1, n
            if (abs(a(r(i), k)) > max_val) then
                max_val = abs(a(r(i), k))
                max_loc = i
            end if
        end do

        call Swap_Integers(r(k), r(max_loc))

        a(r(k+1: n), k) = a(r(k+1: n), k) / a(r(k), k)

        a(r(k+1: n), k+1: n) = a(r(k+1: n), k+1: n) - &
                matmul(a(r(k+1: n), k: k), a(r(k: k), k+1: n))
    end do

    elimination_finish = omp_get_wtime()
end subroutine Fortran9x_Elimination
!!---------------------------------------------------------------------------
subroutine Straightforward_Backsubstition()
    integer :: i, j

    backsubstition_start = omp_get_wtime()

    ! L x = f  =>  x = L \ f
    do i = 1, n
        do j = 1, i-1
            y(r(i)) = y(r(i)) - a(r(i), j) * y(r(j))
        end do
    end do

    ! U y = x  =>  y = U \ x
    do i = n, 1, -1
        do j = i+1, n
            y(r(i)) = y(r(i)) - a(r(i), j) * y(r(j))
        end do

        y(r(i)) = y(r(i)) / a(r(i), i)
    end do

    backsubstition_finish = omp_get_wtime()
end subroutine Straightforward_Backsubstition
!!---------------------------------------------------------------------------
subroutine Fortran9x_Backsubstition()
    integer :: i

    backsubstition_start = omp_get_wtime()

    ! L x = f  =>  x = L \ f
    do i = 1, n
        y(r(i)) = y(r(i)) - dot_product(a(r(i), 1: i-1), y(r(1: i-1)))
    end do

    ! U y = x  =>  y = U \ x
    do i = n, 1, -1
        y(r(i)) = y(r(i)) - dot_product(a(r(i), i+1: n), y(r(i+1: n)))
        y(r(i)) = y(r(i)) / a(r(i), i)
    end do

    backsubstition_finish = omp_get_wtime()
end subroutine Fortran9x_Backsubstition
!!---------------------------------------------------------------------------
end program ge
!!---------------------------------------------------------------------------
