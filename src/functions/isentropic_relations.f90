! Calculates primitive variable vector for given Mach number and stagnation
! conditions
function isentropic_relations(M,V)
  use set_precision, only : prec
  use set_constants, only : one, half
  use set_inputs,    only : gamma, T0, p0, R_gas

  implicit none

  real(prec), dimension(:), intent(in)      :: M
  real(prec), dimension(:,:), intent(out) :: V
  real(prec), dimension(:) :: psi
  real(prec), dimension(:) :: T
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

end function isentropic_relations
