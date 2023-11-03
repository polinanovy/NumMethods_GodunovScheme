program main
use modd

implicit none

integer :: N
real(8) :: C
real(8), allocatable :: u_left(:), u_right(:)
real(8) :: t_stop
real(8) :: g
real(8) :: x_left, x_right

allocate (u_left(3), u_right(3))

! Read the following parameters from file 'INPUT':
call InitializeParameters(u_left, u_right, t_stop, N, C, x_left, x_right)
write(*,*) u_left, u_right, t_stop, N, C, x_left, x_right

end
