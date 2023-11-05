module modd
implicit none
real(8) :: g = 5.d0/3.d0

contains

subroutine InitializeParameters(u_left, u_right, t_stop, N, C, x_left, x_right)
! Subroutine for setting model parameters
integer :: N, i
real(8) :: C, t_stop, x_left, x_right
real(8), allocatable :: u_left(:), u_right(:)
open(unit = 1, file = 'INPUT')
read(1, *) (u_left(i), i=0,2)
read(1, *) (u_right(i), i=0,2)
read(1, *) t_stop
read(1, *) N
read(1, *) C
read(1, *) x_left
read(1, *) x_right
end subroutine InitializeParameters

subroutine CalcLambda(N, rho, v, p, lambda)
!Subroutine for lambda calculation
real(8) :: rho(:), v(:), p(:)
real(8) :: lambda, cs
integer :: N, i
lambda = 0
do i = 0, N-1
 cs = sqrt(g * p(i) / rho(i))
 lambda = max(lambda, abs(v(i)) + cs)
enddo
end subroutine CalcLambda

subroutine InitializeGrid(N, C, x_left, x_right, x, lambda, dx, dt)
! Subroutine for initialising the grid
real(8) :: x_left, x_right, C, lambda
integer :: N, i
real(8) :: dx, dt, x(0:N-1)
dx = (x_right - x_left) / (N - 1)
dt = C * dx / lambda
x(0) = x_left; x(N-1) = x_right
do i = 1, N-2
	x(i) = x(i-1) + dx
enddo
end subroutine InitializeGrid

subroutine CalcEnergy(N, rho, p, E)
!Subroutine for energy calculation
real(8) :: rho(:), p(:)
real(8) :: E(:)
integer :: N, i
do i = 0, N-1
 E(i) = p(i) / (rho(i)*(g-1))
enddo
end subroutine CalcEnergy

subroutine Flux(N, rho, v, p, E, lambda, F_next, F_last)
!Subroutine for flux calculation according to Lax-Friedrichs scheme
integer :: i, N
real(8) :: lambda
real(8) :: rho(:), v(:), p(:), E(:)
real(8) :: F(:,:), U(:,:)
real(8) :: F_next(:,:), F_last(:,:)
!Components of u and F vectors according (6.74)
allocate(F(3,N), U(3,N))
do i = 0, N-1
 F(0,i) = rho(i) * v(i)
 U(0,i) = rho(i)
 F(1,i) = rho(i) * v(i)**2 + p(i)
 U(1,i) = f1(i)
 F(2,i) = rho(i) * v(i) * (E(i) + v(i)**2 / 2.0d0 + p(i) / rho(i))
 U(2,i) = rho(i) * (E(i) + v(i)**2 / 2.0d0)
enddo
!F_next == F_(j+1/2); F_last == F_(j-1/2) in Lax-Friedrichs scheme 
do i = 0, N-2
 F_next(0,i) = 0.5d0 * (F(0,i) + F(0,i+1) - lambda * (U(0,i+1) - U(0,i)))
 F_next(1,i) = 0.5d0 * (F(1,i) + F(1,i+1) - lambda * (U(1,i+1) - U(1,i)))
 F_next(2,i) = 0.5d0 * (F(2,i) + F(2,i+1) - lambda * (U(2,i+1) - U(2,i)))
enddo
F_next(0,N-1) = F(0,N-1)
F_next(1,N-1) = F(1,N-1)
F_next(2,N-1) = F(2,N-1)
F_last(0,0) = F(0,0)
F_last(1,0) = F(1,0)
F_last(2,0) = F(2,0)
do i = 1, N-1
 F_last(0,i) = 0.5d0 * (F(0,i-1) + F(0,i) - lambda * (U(0,i) - U(0,i-1)))
 F_last(1,i) = 0.5d0 * (F(1,i-1) + F(1,i) - lambda * (U(1,i) - U(1,i-1)))
 F_last(2,i) = 0.5d0 * (F(2,i-1) + F(2,i) - lambda * (U(2,i) - U(2,i-1)))
enddo
deallocate(F(3,N), U(3,N))
end subroutine Flux

subroutine GDEquations(N, rho_old, v_old, E_old, dt, dx, F_next, F_last, rho, v, E, p)
!Subroutine for density, velocity and energy calculation according (6.70)-(6.73)
real(8) :: dt, dx
real(8) :: rho_old(:), v_old(:), E_old(:), F_next(:), F_last(:)
real(8) :: rho(:), v(:), E(:), p(:)
integer :: N, i
do i = 0, N-1
 rho(i) = rho_old(i) - dt * (F_next(0, i) - F_last(0, i)) / dx
 v(i) = rho_old(i) * v_old(i) / rho(i) - dt * (F_next(1, i) - F_last(1, i)) / (dx * rho(i))
 E(i) = -v(i)**2 / 2.0d0 + rho_old(i) * (E_old(i) + v_old(i)**2 / 2.0d0) / rho(i) &
       &- dt * (F_next(2, i) - F_last(2, i)) / (dx * rho(i))
 p(i) = rho(i) * E(i) * (g - 1)
enddo
end subroutine GDEquations

subroutine SaveData(N, rho, v, p, x, t_stop, C)
real(8) :: x(0:N-1), rho(0:N-1), v(0:N-1), p(0:N-1), C, t_stop
integer :: N, i
open(unit = 2, file = 'RESULT')
write(2,*) N, C, t_stop
do i = 0, N-1
	write(2,*) x(i), rho(i), v(i), p(i)
enddo
end subroutine SaveData

end module
