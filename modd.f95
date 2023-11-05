module modd
implicit none

contains

subroutine velocities(left, right, g, c_L, c_R, v_A, v_B)
real(8) :: g, left(:), right(:)
real(8) :: c_L, c_R, v_A, v_B
c_L = sqrt(g*left(3) / left(1))
c_R = sqrt(g*right(3) / right(1))
v_A = 2*c_L / (g-1) * ((right(3)/left(3))**((g-1)/(2*g)) - 1)
v_B = 1 / (right(1)*c_R) * (left(3)-right(3)) / sqrt((g+1)/(2*g)*left(3)/right(3)+(g-1)/(2*g))
end subroutine velocities

subroutine F_A(p, c_L, c_R, v_A, v_B, g, left, right, k, v)
real(8) :: c_L, c_R, v_A, v_B, g, left(:), right(:)
real(8) :: v, p, v_L, v_R
integer :: k
if (k==1) then
 v_L = left(2) - 2*c_L / (g-1) * ((p/left(3))**((g-1)/(2*g)) - 1)
 v_R = right(2) + 1 / (right(1)*c_R) * (p-right(3)) / sqrt((g+1)/(2*g)*p/right(3)+(g-1)/(2*g))
else if (k==2) then
 v_L = left(2) - 1/(left(1)*c_L)*(p-left(3))/sqrt((g+1)/(2*g)*p/left(3)+(g-1)/(2*g))
 v_R = right(2) + 1/(right(1)*c_R)*(p-right(3))/sqrt((g+1)/(2*g)*p/right(3)+(g-1)/(2*g))
else if (k==3) then
 v_L = left(2) - 2*c_L / (g-1) * ((p/left(3))**((g-1)/(2*g)) - 1)
 v_R = right(2) + 2*c_R/(g-1)*((p/right(3))**((g-1)/(2*g))-1)
end if
v = v_L - v_R
end subroutine F_A

subroutine newton(left, right, g, c_L, c_R, v_A, v_B, k, p_s)
real(8) :: left(:), right(:), g, c_L, c_R, v_A, v_B
real(8) :: p0, p1, p2, F, dF, temp, F1, F2, p_s
real(8), parameter :: eps = 1e-4, dp = 1e-3 !eps-точность; dp-приращение
integer :: k, i
p0 = 0.5d0 * (left(3)+right(3)) !начальное приближение
p1 = p0
p2 = 0.0
temp = 1.0
i = 0
do while (abs(p2-temp) > eps)
 temp = p1
 call F_A(p1+dp, c_L, c_R, v_A, v_B, g, left, right, k, F1)
 call F_A(p1-dp, c_L, c_R, v_A, v_B, g, left, right, k, F2)
 dF = (F1 - F2) / (2*dp)
 call F_A(p1, c_L, c_R, v_A, v_B, g, left, right, k, F)
 p2 = p1 - F/dF
 p1 = p2
 i = i + 1
end do
p_s = p2
end subroutine newton

subroutine configuration_A(left, right, g, c_L, c_R, v_A, v_B, x, t, RHO, V, P)
real(8) :: left(:), right(:), g, c_L, c_R, v_A, v_B, x(:), t
real(8) :: p_s, v_L_s, rho_L_s, v_R_s, rho_R_s, v_H_RW, v_T_RW, D, RHO(:), V(:), P(:), x_H_RW, x_T_RW, x_s, x_D
real(8), dimension(100) :: v_RW, rho_RW, p_RW
integer :: i, k
k = 1
call newton(left, right, g, c_L, c_R, v_A, v_B, k, p_s)
!Область между ВР и УВ, содержащая КР:
v_L_s = left(2) - 2*c_L / (g-1) * ((p_s/left(3))**((g-1)/(2*g)) - 1)
rho_L_s = left(1) * (p_s/left(3))**(1/g)
v_R_s = right(2) + 1 / (right(1)*c_R) * (p_s-right(3)) / sqrt((g+1)/(2*g)*p_s/right(3)+(g-1)/(2*g))
rho_R_s = right(1) * ((g+1)*p_s/right(3)+(g-1)) / ((g+1)+(g-1)*p_s/right(3))
!rho_R_s = right(1) * ((g+1)*right(3)/p_s+(g-1)) / ((g+1)+(g-1)*right(3)/p_s) !была опечатка
x_s = t * v_L_s
!ВР разрежения, распространяющаяся влево:
v_H_RW = left(2) - c_L !Скорость движения переднего фронта
x_H_RW = t * v_H_RW
v_T_RW = v_R_s - c_L*(p_s/left(3))**((g-1)/(2*g)) !Скорость движения заднего фронта
x_T_RW = t * v_T_RW
do i = 1, 100
 v_RW(i) = left(2) * (g-1)/(g+1) + 2/(g+1)*(x(i)/t+c_L)
 rho_RW(i) = left(1) * (2/(g+1)+(g-1)/((g+1)*c_L)*(left(2)-x(i)/t))**(2/(g-1)) !Скорость, плотность и давление газа внутри ВР
 p_RW(i) = left(3) * (2/(g+1)+(g-1)/((g+1)*c_L)*(left(2)-x(i)/t))**(2*g/(g-1))
end do
!УВ, распространяющаяся вправо:
D = right(2) + c_R*sqrt((g+1)/(2*g)*p_s/right(3)+(g-1)/(2*g))
x_D = t * D
!Все вместе:
do i = 1, 200
 if (x(i)<x_H_RW) then
  RHO(i) = left(1); V(i) = left(2); P(i) = left(3)
 else if (x(i)>x_H_RW .and. x(i)<x_T_RW) then
  RHO(i) = rho_RW(i); V(i) = v_RW(i); P(i) = p_RW(i)
 else if (x(i)>x_T_RW .and. x(i)<x_s) then
  RHO(i) = rho_L_s; V(i) = v_L_s; P(i) = p_s
 else if (x(i)>x_s .and. x(i)<x_D) then
  RHO(i) = rho_R_s; V(i) = v_R_s; P(i) = p_s
 else if (x(i)>x_D) then
  RHO(i) = right(1); V(i) = right(2); P(i) = right(3)
 end if
end do
end subroutine configuration_A

subroutine configuration_B(left, right, g, c_L, c_R, v_A, v_B, x, t, RHO, V, P)
real(8) :: left(:), right(:), g, c_L, c_R, v_A, v_B, x(:), t
real(8) :: p_s, v_L_s, rho_L_s, v_R_s, x_s, rho_R_s, D_L, x_D_L, D_R, x_D_R, RHO(:), V(:), P(:)
integer :: i, k
k = 2
call newton(left, right, g, c_L, c_R, v_A, v_B, k, p_s)
!Область между фронтами УВ:
v_L_s = left(2) - 1/(left(1)*c_L)*(p_s-left(3))/sqrt((g+1)/(2*g)*p_s/left(3)+(g-1)/(2*g))
rho_L_s = left(1) * ((g+1)*p_s/left(3)+(g-1))/((g-1)*p_s/left(3)+(g+1))
v_R_s = right(2) + 1/(right(1)*c_R)*(p_s-right(3))/sqrt((g+1)/(2*g)*p_s/right(3)+(g-1)/(2*g))
rho_R_s = right(1) * ((g+1)*p_s/right(3)+(g-1))/((g-1)*p_s/right(3)+(g+1))
x_s = t * v_L_s
!Скорости фронтов УВ:
D_L = left(2) - c_L*sqrt((g+1)/(2*g)*p_s/left(3)+(g-1)/(2*g))
x_D_L = t * D_L
D_R = right(2) + c_R*sqrt((g+1)/(2*g)*p_s/right(3)+(g-1)/(2*g))
x_D_R = t* D_R
!Все вместе:
do i = 1, 200
 if (x(i)<x_D_L) then
  RHO(i) = left(1); V(i) = left(2); P(i) = left(3)
 else if (x(i)>x_D_L .and. x(i)<x_s) then
  RHO(i) = rho_L_s; V(i) = v_L_s; P(i) = p_s
 else if (x(i)>x_s .and. x(i)<x_D_R) then
  RHO(i) = rho_R_s; V(i) = v_R_s; P(i) = p_s
 else if (x(i)>x_D_R) then
  RHO(i) = right(1); V(i) = right(2); P(i) = right(3)
 end if
end do
end subroutine configuration_B

subroutine configuration_C(left, right, g, c_L, c_R, v_A, v_B, x, t, RHO, V, P)
real(8) :: left(:), right(:), g, c_L, c_R, v_A, v_B, x(:), t
real(8) :: p_s, v_L_s, rho_L_s, v_R_s, rho_R_s, x_s, v_H_RW_L, v_T_RW_L
real(8) :: x_H_RW_L, x_T_RW_L, v_H_RW_R, v_T_RW_R, x_H_RW_R, x_T_RW_R, RHO(:), V(:), P(:)
real(8), dimension(100) :: v_RW_L, rho_RW_L, p_RW_L, v_RW_R, rho_RW_R, p_RW_R
integer :: i, k
k = 3
call newton(left, right, g, c_L, c_R, v_A, v_B, k, p_s)
!Область между ВР и УВ, содержащая КР:
v_L_s = left(2) - 2*c_L / (g-1) * ((p_s/left(3))**((g-1)/(2*g)) - 1)
rho_L_s = left(1) * (p_s/left(3))**(1/g)
v_R_s = right(2) + 2*c_R/(g-1)*((p_s/right(3))**((g-1)/(2*g))-1)
rho_R_s = right(1) * (p_s/right(3))**(1/g)
x_s = t * v_L_s
!ВР разрежения, распространяющаяся влево:
v_H_RW_L = left(2) - c_L !Скорость движения переднего фронта
x_H_RW_L = t * v_H_RW_L
v_T_RW_L = v_R_s - c_L*(p_s/left(3))**((g-1)/(2*g)) !Скорость движения заднего фронта
x_T_RW_L = t * v_T_RW_L
do i = 1, 100
 v_RW_L(i) = left(2) * (g-1)/(g+1) + 2/(g+1)*(x(i)/t+c_L)
 rho_RW_L(i) = left(1) * (2/(g+1)+(g-1)/((g+1)*c_L)*(left(2)-x(i)/t))**(2/(g-1)) !Скорость, плотность и давление газа внутри ВР
 p_RW_L(i) = left(3) * (2/(g+1)+(g-1)/((g+1)*c_L)*(left(2)-x(i)/t))**(2*g/(g-1))
end do
!ВР разрежения, распространяющаяся вправо:
v_H_RW_R = right(2) + c_R
x_H_RW_R = t * v_H_RW_R
v_T_RW_R = v_R_s + c_R*(p_s/right(3))**((g-1)/(2*g))
x_T_RW_R = t * v_T_RW_R
do i = 1, 100
 v_RW_R(i) = right(2) * (g-1)/(g+1) + 2/(g+1)*(x(i+100)/t-c_R)
 rho_RW_R(i) = right(1) * (2/(g+1)-(g-1)/((g+1)*c_R)*(right(2)-x(i+100)/t))**(2/(g-1)) !Скорость, плотность и давление газа внутри ВР
 p_RW_R(i) = right(3) * (2/(g+1)-(g-1)/((g+1)*c_R)*(right(2)-x(i+100)/t))**(2*g/(g-1))
end do
!Все вместе:
do i = 1, 200
 if (x(i)<x_H_RW_L) then
  RHO(i) = left(1); V(i) = left(2); P(i) = left(3)
 else if (x(i)==x_H_RW_L) then
  V(i) = v_H_RW_L; RHO(i) = left(1); P(i) = left(3)
 else if (x(i)>x_H_RW_L .and. x(i)<x_T_RW_L) then
  RHO(i) = rho_RW_L(i); V(i) = v_RW_L(i); P(i) = p_RW_L(i)
 !else if (x(i)==x_T_RW_L) then
  !V(i) = v_T_RW_L; RHO(i) = rho_L_s; P(i) = p_s
 else if (x(i)>x_T_RW_L .and. x(i)<x_T_RW_R) then
  RHO(i) = rho_L_s; V(i) = v_L_s; P(i) = p_s
 !else if (x(i)==x_T_RW_R) then
  !V(i) = v_T_RW_R; RHO(i) = rho_R_s; P(i) = p_s
 else if (x(i)>x_T_RW_R .and. x(i)<x_H_RW_R) then
  RHO(i) = rho_RW_R(i-100); V(i) = v_RW_R(i-100); P(i) = p_RW_R(i-100)
 !else if (x(i)==x_H_RW_R) then
  !V(i) = v_H_RW_R; RHO(i) = right(1); P(i) = right(3)
 else if (x(i)>x_H_RW_R) then
  RHO(i) = right(1); V(i) = right(2); P(i) = right(3)
 end if
end do
end subroutine configuration_C

end module
