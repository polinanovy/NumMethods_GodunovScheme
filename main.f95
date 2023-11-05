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

allocate (u_left(3), u_right(3))

! Read the following parameters from file 'INPUT':
call InitializeParameters(u_left, u_right, t_stop, N, C, x_left, x_right)

allocate(x(N))
! Get grid step sizes dx and dt, as well as the array of coordinates "x":
call InitializeGrid(N, C, x_left, x_right, x, lambda, dx, dt) 

allocate(rho(N), v(N), p(N), E(N))
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

! Start timer:
t = dt

allocate(F_next(3,N), F_last(3,N), rho_new(N), v_new(N), E_new(N), p_new(N))
do while (t <= t_stop)
 call CalcLambda(N, rho, v, p, lambda)
 call Flux(N, rho, v, p, E, lambda, F_next, F_last)
 call GDEquations(N, rho, v, E, dt, dx, F_next, F_last, rho_new, v_new, E_new, p_new)
 do i = 0, N
  rho(i) = rho_new(i)
  v(i) = v_new(i)
  E(i) = E_new(i)
  p(i) = p_new(i)
 enddo
 ! Update timer:
 t = t + dt
end do

! Save data to file 'RESULT':
call SaveData(N, rho, v, p, x, t_stop, C)

deallocate(u_left, u_right, x, rho, v, p, E, F_next, F_last, rho_new, v_new, E_new, p_new)

end
