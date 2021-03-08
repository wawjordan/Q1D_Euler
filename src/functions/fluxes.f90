module fluxes
  
  use set_precision, only : prec
  use set_constants, only : one, three, half
  use fluid_constants, only : gamma
  
  implicit none
  
contains
  
subroutine central_flux(U,F)
  
  real(prec), dimension(:,:), intent(inout)  :: U
  real(prec), dimension(:,:), intent(inout)  :: F
  
  F(:,1) = U(:,2)
  F(:,2) = half*(three-gamma)*( U(:,2)**2 )/U(:,1) + (gamma-one)*U(:,3)
  F(:,3) = U(:,3)*U(:,2)/U(:,1) + U(:,2)/U(:,1)*( (gamma-one)*U(:,3) - half*(gamma-one)*U(:,2)**2/U(:,1) )
  
end subroutine central_flux

end module fluxes
