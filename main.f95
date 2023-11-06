program main
use modd

implicit none

integer :: N, i
real(8) :: C
real(8), allocatable :: u_left(:), u_right(:)
real(8) :: dx, dt, t, t_stop
real(8) :: x_left, x_right
real(8) :: lambda
real(8), allocatable :: rho(:), v(:), p(:), E(:), x(:)
real(8), allocatable :: F_next(:,:), F_last(:,:), rho_new(:), v_new(:), E_new(:), p_new(:)

allocate(u_left(0:2), u_right(0:2))

! Read the following parameters from file 'INPUT':
call InitializeParameters(u_left, u_right, t_stop, N, C, x_left, x_right)

allocate(rho(0:N-1), v(0:N-1), p(0:N-1), E(0:N-1))
do i = 0, N/2-1
	rho(i) = u_left(0)
	v(i)   = u_left(1)
	p(i)   = u_left(2)
enddo
do i = N/2, N-1
	rho(i) = u_right(0)
	v(i)   = u_right(1)
	p(i)   = u_right(2)
enddo
call CalcEnergy(N, rho, p, E)

! Get grid step sizes dx and dt, as well as the array of coordinates "x":
allocate(x(0:N-1))
call CalcLambda(N, rho, v, p, lambda)
call InitializeGrid(N, C, x_left, x_right, x, lambda, dx, dt) 

! Start timer:
t = 0.d0

allocate(F_next(0:2,0:N-1), F_last(0:2,0:N-1), rho_new(0:N-1), v_new(0:N-1), E_new(0:N-1), p_new(0:N-1))
do while (t <= t_stop)
	call CalcLambda(N, rho, v, p, lambda)
	call TimeStep(N, C, lambda, dx, dt)
	call Flux(N, rho, v, p, E, lambda, F_next, F_last)
	call GDEquations(N, rho, v, E, dt, dx, F_next, F_last, rho_new, v_new, E_new, p_new, u_left, u_right)
	rho = rho_new
	v = v_new
	E = E_new
	p = p_new
	! Update timer:
	t = t + dt
end do

! Save data to file 'RESULT':
call SaveData(N, rho, v, p, x, t_stop, C)

deallocate(u_left, u_right, x, rho, v, p, E, F_next, F_last, rho_new, v_new, E_new, p_new)

end
