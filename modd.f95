module modd
use matrix_solver
implicit none
real(8) :: pi = 4 * atan(1.d0)

contains

subroutine InitializeParameters(u_left, u_right, t_stop, N, C)
! Subroutine for setting model parameters
integer :: N, i
real(8) :: C, t_stop
real(8), allocatable :: u_left(:), u_right(:)
open(unit = 1, file = 'INPUT')
read(1, *) (u_left(i), i=1,3)
read(1, *) (u_right(i), i=1,3)
read(1, *) t_stop
read(1, *) N
read(1, *) C
end subroutine InitializeParameters

end module
