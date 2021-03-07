module variable_conversion
   
  use set_precision,   only : prec
  use set_constants,   only : one, half
  use fluid_constants, only : gamma, R_gas
  use set_inputs,      only : p0,T0

  implicit none

contains

subroutine speed_of_sound(pressure,rho,sound_speed)

  real(prec), dimension(:), intent(in)  :: pressure, rho
  real(prec), dimension(:), intent(out) :: sound_speed
  
  sound_speed = sqrt(gamma*pressure/rho)

end subroutine speed_of_sound


subroutine prim2cons(U,V)

  real(prec), dimension(:,:), intent(inout)  :: U
  real(prec), dimension(:,:), intent(in)     :: V

  U(:,1) = V(:,1)
  U(:,2) = V(:,1)*V(:,2)
  U(:,3) = ( V(:,3)/(gamma - one) ) + half*V(:,1)*V(:,2)**2

end subroutine prim2cons


subroutine cons2prim(U,V)

  real(prec), dimension(:,:), intent(in)    :: U
  real(prec), dimension(:,:), intent(inout) :: V

  V(:,1) = U(:,1)
  V(:,2) = U(:,2)/U(:,1)
  V(:,3) = (gamma - one)*U(:,3) - half*(gamma - one)*U(:,2)**2/U(:,1)

end subroutine cons2prim


subroutine isentropic_relations(M,V,T)

  real(prec), dimension(:,:), intent(inout) :: V
  real(prec), dimension(:),   intent(in) :: M
  real(prec), dimension(:),   intent(out) :: T
  
  T(:) = T0/(one + half*(gamma - one)*M(:)**2)
  V(:,3) = 1000.0_prec*p0/(one + half*(gamma - one)*M(:)**2)**(gamma/(gamma-1))
  V(:,1) = V(:,3)/(R_gas*T(:))
  call speed_of_sound(V(:,3),V(:,1),V(:,2))
  V(:,2) = M(:)*V(:,2)

end subroutine isentropic_relations


end module variable_conversion
