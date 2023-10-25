program main
use modd

implicit none

integer :: N
real(8) :: C
real(8), allocatable :: u_left(:), u_right(:)
real(8) :: t_stop
real(8) :: g

g = 5.d0/3.d0

! Read the following parameters from file 'INPUT':
call InitializeParameters(u_left, u_right, t_stop, N, C)

end
