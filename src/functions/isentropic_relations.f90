! Calculates primitive variable vector for given Mach number and stagnation
! conditions
subroutine isentropic_relations(M,V,T)

  use set_precision,   only : prec
  use set_constants,   only : one, half
  use fluid_constants, only : gamma, R_gas
  use set_inputs,      only : T0, p0
  use speed_of_sound

  implicit none

  real(prec), dimension(:,:), intent(inout) :: V
  real(prec), dimension(:),   intent(inout) :: M
  real(prec), dimension(:),   intent(inout) :: T

  real(prec), dimension(:) :: psi
  real(prec), dimension(:) :: p
  real(prec), dimension(:) :: rho
  real(prec), dimension(:) :: a
  real(prec), dimension(:) :: u
  
  psi = one + half*(gamma - one)*M**2
  T   = T0/psi
  p   = p0/psi**(gamma/(gamma - one))
  rho = p/(R_gas*T)
  a   = speed_of_sound(p,rho)
  u   = M*a

  V(:,1) = rho
  V(:,2) = u
  V(:,3) = p

end subroutine isentropic_relations
