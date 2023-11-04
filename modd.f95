module modd
implicit none
real(8) :: pi = 4 * atan(1.d0)
real(8) :: g = 5.d0/3.d0

contains

subroutine InitializeParameters(u_left, u_right, t_stop, N, C, x_left, x_right)
! Subroutine for setting model parameters
integer :: N, i
real(8) :: C, t_stop, x_left, x_right
real(8), allocatable :: u_left(:), u_right(:)
open(unit = 1, file = 'INPUT')
read(1, *) (u_left(i), i=1,3)
read(1, *) (u_right(i), i=1,3)
read(1, *) t_stop
read(1, *) N
read(1, *) C
read(1, *) x_left
read(1, *) x_right
end subroutine InitializeParameters

subroutine InitializeGrid(N, C, D, x_left, x_right, x, lambda, dx, dt)   !непонятно, что с t делать
! Subroutine for initialising the grid
real(8) :: x_left, x_right, C, D
integer :: N, i
real(8) :: dx, dt, x(0:N-1)
dx = (x_right - x_left) / (N - 1)
dt = C * dx / lambda
x(0) = x_left; x(N-1) = x_right
do i = 1, N-2
	x(i) = x(i-1) + dx
enddo
end subroutine InitializeGrid

subroutine Scheme()

end subroutine Scheme

subroutine CalcLambda(N, rho, v, p, lambda)
!Subroutine for lambda calculation
real(8) :: rho(:), v(:), p(:)
real(8) :: lambda, cs
integer :: N, i
lambda = 0
do i = 0, N-1
 cs = sqrt(g * p(i) / rho(i))
 max(lambda, abs(v(i)) + cs)
enddo
end subroutine CalcLambda

end module
